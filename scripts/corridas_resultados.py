# Corridas del capítulo 7 con el núcleo del simulador v1 (sin servidor ni pantalla).
# Tres tipos: referencia, perturbación y cable ruidoso; 5 repeticiones cada uno, con semillas fijas.
import importlib.util, os, tempfile, random, statistics as st, json, math

ESC = os.path.join(tempfile.gettempdir(), "pendulo_corridas") + os.sep  # carpeta temporal para las bases de prueba
os.makedirs(ESC, exist_ok=True)
spec = importlib.util.spec_from_file_location("simulador", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulador.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
REPS = 5

def nueva(semilla, ruidoso=False):
    random.seed(semilla)
    db = ESC + "tmp.db"
    for f in (db, db + "-wal", db + "-shm"):
        if os.path.exists(f): os.remove(f)
    s = m.Simulador(db); s._abrir_db(); s._nueva_sesion(); s.cable_ruidoso = ruidoso
    s.desvios = []                                      # desvío torre-plano en cada medición con señal
    orig = s.arduino._medir_y_corregir
    def espia(t, ms, az):                               # registra el desvío real antes de cada medición
        if s.arduino.estado == "SIGUIENDO":
            s.desvios.append(m.desvio(az, s.arduino.az_torre))
        orig(t, ms, az)
    s.arduino._medir_y_corregir = espia
    return s

def hasta(s, horas):
    while s.t < horas * 3600:
        s._avanzar(min(600, horas * 3600 - s.t + 1e-6))

def periodo(s):
    i = s.ingestor
    return 2.0 * (i.ms_ultimo_cruce - i.ms_primer_cruce) / 1000.0 / (i.cruces - 1)

R = {}
# 1) Referencia: 24 h sin fallas
ref = []
for k in range(REPS):
    s = nueva(100 + k); hasta(s, 24)
    om, se = s.ingestor.omega()
    ref.append(dict(periodo=periodo(s), omega=om, u=se, desvio_max=max(abs(d) for d in s.desvios), huecos=s.ingestor.huecos))
    s.con.close()
R["referencia"] = ref

# 2) Perturbación: 2 h normales, giro de 25°, luego 6 h más
per = []
for k in range(REPS):
    s = nueva(200 + k); hasta(s, 2)
    s._perturbar(); t_p = s.t; n_antes = len(s.desvios)
    while s.ingestor.huecos == 0: s._avanzar(5.0)        # hasta que se cierra el hueco (recuperación)
    t_r = s.t
    row = s.con.execute("SELECT inicio_ms, fin_ms, pasos_barrido FROM huecos").fetchone()
    i_rec = len(s.desvios)
    hasta(s, t_r / 3600 + 6)
    om, se = s.ingestor.omega()
    despues = [abs(d) for d in s.desvios[i_rec + 3:]]    # tras los primeros ciclos del tramo nuevo
    per.append(dict(seg_hasta_perdida=(row[0] / 1000 - t_p), hueco_s=(row[1] - row[0]) / 1000, pasos_barrido=row[2],
                    seg_total=t_r - t_p, desvio_max_despues=max(despues), omega_tramo2=om, u_tramo2=se, tramo=s.ingestor.tramo_numero))
    s.con.close()
R["perturbacion"] = per

# 3) Cable ruidoso (2 % de líneas dañadas) contra la misma semilla sin ruido, 24 h
cab = []
for k in range(REPS):
    s0 = nueva(300 + k, False); hasta(s0, 24); om0, se0 = s0.ingestor.omega(); s0.con.close()
    s1 = nueva(300 + k, True); hasta(s1, 24); om1, se1 = s1.ingestor.omega()
    i = s1.ingestor
    cab.append(dict(omega_sin=om0, omega_con=om1, dif=abs(om1 - om0), descartadas=i.lineas_descartadas,
                    total=i.lineas_ok + i.lineas_descartadas, incompletos=i.ciclos_incompletos, huecos=i.huecos))
    s1.con.close()
R["cable"] = cab
json.dump(R, open(ESC + "resultados.json", "w"), indent=1)

f = lambda xs: (st.mean(xs), st.stdev(xs), min(xs), max(xs))
print("Ω teórica", m.OMEGA_TEORICA, "período teórico", 2 * math.pi * math.sqrt(m.LARGO / 9.81))
for k, v in R.items():
    print("==", k)
    for c in v[0]:
        xs = [r[c] for r in v]
        if isinstance(xs[0], (int, float)) and xs[0] is not None:
            print(f"  {c}: media {f(xs)[0]:.6g} desv {f(xs)[1]:.3g} min {f(xs)[2]:.6g} max {f(xs)[3]:.6g}")
