# Corridas con ruido de 2,7° por lectura (RUIDO_PICO = 0,4 V): referencia 24 h y perturbación, 5 repeticiones cada una.
import importlib.util, os, tempfile, random, statistics as st, math, json
ESC = os.path.join(tempfile.gettempdir(), "pendulo_corridas") + os.sep  # carpeta temporal para las bases de prueba
os.makedirs(ESC, exist_ok=True)
spec = importlib.util.spec_from_file_location("simulador", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulador.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.RUIDO_PICO = 0.4
print("grados por lectura:", round(math.sqrt(2) * m.RUIDO_PICO / m.PENDIENTE, 3))
def nueva(semilla):
    random.seed(semilla); db = ESC + "tmp27.db"
    for f in (db, db + "-wal", db + "-shm"):
        if os.path.exists(f): os.remove(f)
    s = m.Simulador(db); s._abrir_db(); s._nueva_sesion(); s.desvios = []
    orig = s.arduino._medir_y_corregir
    def espia(t, ms, az):
        if s.arduino.estado == "SIGUIENDO": s.desvios.append(m.desvio(az, s.arduino.az_torre))
        orig(t, ms, az)
    s.arduino._medir_y_corregir = espia
    return s
def hasta(s, h):
    while s.t < h * 3600: s._avanzar(min(600, h * 3600 - s.t + 1e-6))
ref = []
for k in range(5):
    s = nueva(400 + k); hasta(s, 24); om, se = s.ingestor.omega()
    ref.append(dict(omega=om, u=se, desvio_max=max(abs(d) for d in s.desvios), huecos=s.ingestor.huecos)); s.con.close()
per = []
for k in range(5):
    s = nueva(500 + k); hasta(s, 2); s._perturbar(); t_p = s.t
    while s.ingestor.huecos == 0 and s.t < t_p + 3600: s._avanzar(5.0)
    rec = s.ingestor.huecos > 0
    row = s.con.execute("SELECT inicio_ms, fin_ms, pasos_barrido FROM huecos ORDER BY id LIMIT 1").fetchone()
    i_rec = len(s.desvios); t_r = s.t; h0 = s.ingestor.huecos
    hasta(s, t_r / 3600 + 6); om, se = s.ingestor.omega()
    per.append(dict(recupero=rec, seg_perdida=row[0] / 1000 - t_p, hueco_s=(row[1] - row[0]) / 1000, pasos=row[2],
                    desvio_max_despues=max(abs(d) for d in s.desvios[i_rec + 3:]), huecos_extra=s.ingestor.huecos - h0, omega=om, u=se))
    s.con.close()
json.dump(dict(ref=ref, per=per), open(ESC + "resultados_27.json", "w"), indent=1)
for n, v in (("ref", ref), ("per", per)):
    print("==", n)
    for c in v[0]:
        xs = [r[c] for r in v]; print(" ", c, [round(x, 4) if isinstance(x, float) else x for x in xs])
