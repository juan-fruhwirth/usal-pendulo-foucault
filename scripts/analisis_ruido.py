# Análisis de sensibilidad al ruido de los sensores, con el núcleo del simulador v1 (sin servidor ni pantalla).
# Para cada nivel de ruido y cada repetición se simula hasta 24 h y, en varios instantes, se anota:
#   - u(Ω): error estándar de la pendiente que informa la regresión (lo que mostraría el sistema)
#   - Ω medida, para calcular después la dispersión REAL entre repeticiones (error efectivo)
import importlib.util, math, os, random, statistics, json, tempfile

ESC = os.path.join(tempfile.gettempdir(), "pendulo_corridas") + os.sep  # carpeta temporal para las bases de prueba
os.makedirs(ESC, exist_ok=True)
spec = importlib.util.spec_from_file_location("simulador", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulador.py"))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

NIVELES_V = [0.005, 0.05, 0.15, 0.5, 1.0]        # desvío del ruido de cada pico [V]; 0,005 es el valor del v1
HORAS = [0.25, 0.5, 1, 2, 4, 8, 12, 24]            # instantes en que se anota el resultado
REPETICIONES = 5
OBJETIVO = 0.43                                    # °/h

def ruido_en_grados(sigma_pico):
    # V_OE = suma de dos picos -> desvío √2·σ; se pasa a grados con la sensibilidad cerca del centro (V/°)
    return math.sqrt(2) * sigma_pico / m.PENDIENTE

resultados = {}
for sigma in NIVELES_V:
    m.RUIDO_PICO = sigma                            # el simulador lee esta constante en cada medición
    filas = {h: [] for h in HORAS}
    huecos = []
    for rep in range(REPETICIONES):
        random.seed(1000 * rep + int(sigma * 1e4))  # semilla distinta por repetición, reproducible
        db = ESC + "tmp.db"
        for f in (db, db + "-wal", db + "-shm"):
            if os.path.exists(f): os.remove(f)
        s = m.Simulador(db); s._abrir_db(); s._nueva_sesion()
        t_prev = 0.0
        for h in HORAS:
            while s.t < h * 3600:                    # avanza en bloques de 10 min simulados
                s._avanzar(min(600, h * 3600 - s.t + 1e-6))
            om, se = s.ingestor.omega()              # Ω y u(Ω) del tramo actual (sin mínimo de duración)
            filas[h].append((om, se, s.ingestor.tramo_numero))
        huecos.append(s.ingestor.huecos)
        s.con.close()
    resultados[sigma] = (filas, huecos)

teo = m.OMEGA_TEORICA
print(f"Ω teórica = {teo:.4f} °/h | objetivo u(Ω) ≤ {OBJETIVO} °/h | {REPETICIONES} repeticiones por nivel\n")
resumen = []
for sigma, (filas, huecos) in resultados.items():
    gr = ruido_en_grados(sigma)
    print(f"=== Ruido por pico {sigma} V  ≈ {gr:.3f}° por lectura | huecos por corrida: {huecos}")
    print(f"{'horas':>6} | {'u(Ω) informada (media)':>22} | {'dispersión real Ω (desvío)':>26} | {'error medio vs teórica':>22}")
    h_u = h_real = None
    for h in HORAS:
        vals = [f for f in filas[h] if f[0] is not None]
        if len(vals) < 2:
            print(f"{h:>6} | sin datos suficientes"); continue
        u_med = statistics.mean(v[1] for v in vals)
        disp = statistics.stdev(v[0] for v in vals)
        err = statistics.mean(abs(v[0] - teo) for v in vals)
        print(f"{h:>6} | {u_med:>22.4f} | {disp:>26.4f} | {err:>22.4f}")
        if h_u is None and u_med <= OBJETIVO: h_u = h
        if h_real is None and disp <= OBJETIVO and err <= OBJETIVO: h_real = h
    print(f"   -> u(Ω) informada ≤ {OBJETIVO}: desde {h_u} h | error real (dispersión y diferencia) ≤ {OBJETIVO}: desde {h_real} h\n")
    resumen.append(dict(sigma_V=sigma, grados=round(gr, 3), horas_u=h_u, horas_real=h_real, huecos=huecos))
json.dump(resumen, open(ESC + "resumen.json", "w"), indent=1)
