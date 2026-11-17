"""
Simulador del péndulo de Foucault (Pilar, 34,45° S) - versión 1.

Prototipo de gemelo digital del sistema de instrumentación descripto en el Capítulo 4.
Tiene tres partes, que en el sistema real serían tres cosas físicas distintas:

  1. EL MUNDO        - la física del péndulo (con Coriolis). La precesión sale sola.
  2. ARDUINO VIRTUAL - modelo de comportamiento de la bobina de detección, el par de sensores Hall
                       y el motor paso a paso, con la lógica de control de 4.4 (acumulador, pasos de
                       0,18°, pérdida de referencia y barrido). Solo "ve" lo que verían los sensores.
                       Se comunica mandando líneas de texto, como el Arduino real por el cable serie.
  3. INGESTOR        - lo que correría en la Raspberry Pi: lee las líneas, verifica que no estén
                       dañadas, arma cada ciclo y lo guarda en SQLite con el esquema de 4.6
                       (sesiones, tramos, datos_crudos, huecos). Calcula Ω por tramo.

Modelo de comportamiento: los sensores no se simulan a nivel eléctrico (voltajes de bobina,
campo del imán, conversor), sino por su respuesta: cuánto se corrió el plano respecto de la torre,
con ruido y con un rango máximo. El modelo eléctrico detallado queda como línea futura.

Uso:  python -m pip install flask   ->   python simulador.py   ->   abrir http://localhost:5000
"""
import json
import math
import os
import random
import sqlite3
import threading
import time
from collections import deque
from datetime import datetime, timedelta

from flask import Flask, jsonify, request, send_from_directory

# ==============================================================================================
# PARÁMETROS
# ==============================================================================================
# --- Péndulo (diseño) ---
G = 9.80665                      # gravedad [m/s^2]
LARGO = 7.0                      # longitud del cable [m]
LATITUD_GRADOS = -34.45          # Pilar; negativa = hemisferio sur
OMEGA_TIERRA = 7.2921159e-5      # rotación sidérea de la Tierra [rad/s]
AMPLITUD = 0.5                   # amplitud [m] (a definir en el diseño)
DT = 0.02                        # paso de cálculo [s simulados]: detalle interno, no es un sensor
VELOCIDADES = (1, 10, 60, 360, 1000)

# --- Motor y torre (diseño, 4.4) ---
GRADOS_POR_PASO = 0.18           # NEMA 17 + reductora 10:1: 0,18° por paso a la salida
K_M = 1.0 / GRADOS_POR_PASO      # pasos de motor por grado de azimut (~5,56)
D_S = 0.7 * AMPLITUD             # distancia de la torre al centro: ~70 % de la amplitud (antecedente)

# --- Sensores Hall: modelo de comportamiento (SUPUESTOS, a justificar o variar) ---
RANGO_SENSORES = 15.0            # [°] más allá, ningún sensor ve el imán (antecedente, 4.4)
V_MAX = 1.0                      # [V] máximo de V_OE dentro del rango
PENDIENTE = V_MAX * math.pi / RANGO_SENSORES   # [V/°] sensibilidad cerca del centro
K_H = 1.0 / PENDIENTE            # [°/V] calibración ideal. Así incremento = K_H·K_M·V_OE da pasos
RUIDO_PICO = 0.005               # [V] ruido de cada pico
LINEA_BASE = 2.5                 # [V] reposo de la cadena (nivel corrido al medio del conversor)
DERIVA_BASE = 0.02               # [V] deriva térmica de la línea base a lo largo del día
RUIDO_BASE = 0.002               # [V]
S_MAX = 2.0                      # [V] salida del amplificador sumador con el imán centrado
RUIDO_SUMA = 0.02                # [V]
UMBRAL_SUMA = 0.4                # [V] por debajo, los sensores "no ven" el imán
CICLOS_PARA_PERDIDA = 3          # ciclos seguidos sin señal para declarar pérdida de referencia
PASOS_BARRIDO = 50               # pasos por ciclo durante el barrido (9°, menor que la ventana de los sensores)

# --- Reloj, tiempos y cable (SUPUESTOS) ---
RUIDO_RELOJ_MS = 0.5             # [ms] error al sellar el instante del cruce
RETARDO_MEDICION_MS = 2100       # [ms] tras el cruce: el imán ya pasó dos veces por la torre
RETARDO_CORRECCION_MS = 2200     # [ms] tras el cruce: el motor corrige antes del próximo cruce
PROB_CORRUPCION = 0.02           # probabilidad de dañar una línea con "cable con ruido" activado
PERTURBACION_GRADOS = 25.0       # giro brusco del plano con el botón "Perturbar"
MIN_HORAS_OMEGA = 0.25           # duración mínima de un tramo para mostrar Ω en pantalla (15 min)

# --- Derivados ---
W02 = G / LARGO
WZ = OMEGA_TIERRA * math.sin(math.radians(LATITUD_GRADOS))
DWZ = 2.0 * WZ
OMEGA_TEORICA = abs(WZ) * 3600.0 * 180.0 / math.pi   # °/h (≈ 8,509 en Pilar)

BASE = os.path.dirname(os.path.abspath(__file__))
RUTA_DB = os.path.join(BASE, "pendulo.db")
VERSION_ESQUEMA = 1              # se incrementa cada vez que cambian las tablas


# ==============================================================================================
# 1. EL MUNDO: física del péndulo
# ==============================================================================================
def paso_rk4(x, y, vx, vy):
    """Avanza el péndulo un paso DT (Runge-Kutta de 4.º orden).
    x'' = -(g/L) x + 2 Ωz y'      y'' = -(g/L) y - 2 Ωz x'   (el término de Coriolis gira el plano)"""
    a1x = -W02 * x + DWZ * vy
    a1y = -W02 * y - DWZ * vx
    x2, y2 = x + 0.5 * DT * vx, y + 0.5 * DT * vy
    vx2, vy2 = vx + 0.5 * DT * a1x, vy + 0.5 * DT * a1y
    a2x = -W02 * x2 + DWZ * vy2
    a2y = -W02 * y2 - DWZ * vx2
    x3, y3 = x + 0.5 * DT * vx2, y + 0.5 * DT * vy2
    vx3, vy3 = vx + 0.5 * DT * a2x, vy + 0.5 * DT * a2y
    a3x = -W02 * x3 + DWZ * vy3
    a3y = -W02 * y3 - DWZ * vx3
    x4, y4 = x + DT * vx3, y + DT * vy3
    vx4, vy4 = vx + DT * a3x, vy + DT * a3y
    a4x = -W02 * x4 + DWZ * vy4
    a4y = -W02 * y4 - DWZ * vx4
    c = DT / 6.0
    return (x + c * (vx + 2 * vx2 + 2 * vx3 + vx4),
            y + c * (vy + 2 * vy2 + 2 * vy3 + vy4),
            vx + c * (a1x + 2 * a2x + 2 * a3x + a4x),
            vy + c * (a1y + 2 * a2y + 2 * a3y + a4y))


def desvio(az_plano, az_torre):
    """Diferencia entre el plano y la torre, en (-90°, 90°]. El plano es una recta: 0° y 180° son lo mismo."""
    return (az_plano - az_torre + 90.0) % 180.0 - 90.0


# ==============================================================================================
# 2. ARDUINO VIRTUAL: modelo de comportamiento de sensores y motor + lógica de control (4.4)
# ==============================================================================================
class ArduinoVirtual:
    """Recibe del mundo solo lo que el hardware podría percibir en cada cruce por el centro, y
    responde con líneas de texto, igual que el firmware real por el puerto serie."""

    def __init__(self, enviar, programar):
        self.enviar = enviar             # función que pone una línea en el "cable"
        self.programar = programar       # agenda una tarea para un instante futuro (temporizador del Arduino)
        self.az_torre = 0.0              # posición física de la torre [°]: la mueve el motor
        self.acumulador = 0.0            # fracción de paso pendiente (4.4, "Acumulador de corrección")
        self.pasos_acumulados = 0        # pasos ejecutados en el tramo actual
        self.estado = "SIGUIENDO"        # SIGUIENDO | BARRIDO
        self.ciclos_sin_senal = 0
        self.pasos_barrido = 0

    def al_detectar_cruce(self, t, az_plano_real, rumbo):
        """Bobina de detección: el imán cruzó el centro. az_plano_real y rumbo son la 'verdad' del mundo;
        el Arduino solo los usa a través de lo que medirían los sensores."""
        ms = round(t * 1000 + random.gauss(0, RUIDO_RELOJ_MS))     # instante sellado por el Arduino [ms]
        # ¿En este medio ciclo la esfera va hacia el lado de la torre? Solo entonces la sobrevuela.
        hacia_torre = math.cos(math.radians(rumbo - self.az_torre)) > 0
        self.enviar(f"C;{ms};{1 if hacia_torre else -1}")
        if not hacia_torre:
            return                                                   # este medio ciclo no se mide
        # La medición termina cuando el imán ya pasó dos veces por la torre: se agenda para ese momento,
        # así las líneas salen por el cable a su hora y no antes.
        self.programar(t + RETARDO_CORRECCION_MS / 1000.0,
                       lambda: self._medir_y_corregir(t, ms, az_plano_real))

    def _medir_y_corregir(self, t, ms, az_plano_real):
        """Ventana de muestreo + corrección (4.4), ejecutada al cierre de la ventana."""
        # --- Sensores Hall (modelo de comportamiento) ---
        d = desvio(az_plano_real, self.az_torre)                    # cuánto se corrió el plano [°]
        visible = abs(d) < RANGO_SENSORES
        v_oe_real = V_MAX * math.sin(math.pi * d / RANGO_SENSORES) if visible else 0.0   # crece y decae a 0 en ±15°
        suma = (S_MAX * math.cos(math.pi * d / (2 * RANGO_SENSORES)) if visible else 0.0) + random.gauss(0, RUIDO_SUMA)
        linea_base = LINEA_BASE + DERIVA_BASE * math.sin(2 * math.pi * t / 86400) + random.gauss(0, RUIDO_BASE)
        pico_ida = v_oe_real / 2 + random.gauss(0, RUIDO_PICO)      # ya referido a la línea base
        pico_regreso = v_oe_real / 2 + random.gauss(0, RUIDO_PICO)  # misma recta: mismo signo que la ida
        v_oe = pico_ida + pico_regreso
        v_d = pico_ida - pico_regreso
        ms_m = ms + RETARDO_MEDICION_MS
        ms_k = ms + RETARDO_CORRECCION_MS

        # --- Barrido: buscando de nuevo al imán ---
        if self.estado == "BARRIDO":
            if suma < UMBRAL_SUMA:
                self.az_torre += PASOS_BARRIDO * GRADOS_POR_PASO     # sigue girando la torre
                self.pasos_barrido += PASOS_BARRIDO
                return                                               # durante el hueco no hay datos
            # lo reencontró: nuevo tramo con origen en cero (4.4, "Recuperación de referencia")
            self.enviar(f"R;{ms};{self.pasos_barrido}")
            self.estado = "SIGUIENDO"
            self.acumulador, self.pasos_acumulados = 0.0, 0
            self.pasos_barrido, self.ciclos_sin_senal = 0, 0

        # --- Medición ---
        self.enviar(f"M;{ms_m};{linea_base:.4f};{pico_ida:.4f};{pico_regreso:.4f};{v_oe:.4f};{v_d:.4f};{suma:.4f}")

        # --- Sin señal: acumulador congelado; tras varios ciclos, pérdida de referencia ---
        if suma < UMBRAL_SUMA:
            self.ciclos_sin_senal += 1
            self.enviar(f"K;{ms_k};0.0000;0;{self.pasos_acumulados};{self.acumulador:.4f}")
            if self.ciclos_sin_senal >= CICLOS_PARA_PERDIDA:
                self.enviar(f"P;{ms_k + 1};{self.pasos_acumulados};{self.acumulador:.4f}")
                self.estado = "BARRIDO"
            return
        self.ciclos_sin_senal = 0

        # --- Corrección (4.4): incremento = K_H · K_M · V_OE, solo pasos enteros ---
        incremento = K_H * K_M * v_oe
        self.acumulador += incremento
        pasos = int(self.acumulador)                 # parte entera con signo: el motor no da medio paso
        if pasos != 0:
            self.az_torre += pasos * GRADOS_POR_PASO # el motor gira la torre
            self.acumulador -= pasos
            self.pasos_acumulados += pasos
        self.enviar(f"K;{ms_k};{incremento:.4f};{pasos};{self.pasos_acumulados};{self.acumulador:.4f}")


# ==============================================================================================
# 3. INGESTOR: lo que correría en la Raspberry Pi (lee el cable, valida, guarda, analiza)
# ==============================================================================================
CAMPOS_POR_TIPO = {"C": 3, "M": 8, "K": 6, "P": 4, "R": 3}
NOMBRES = {"C": "Cruce", "M": "Medición", "K": "Corrección", "P": "Pérdida de referencia", "R": "Recuperación"}


def checksum(texto):
    """XOR de todos los caracteres: detecta líneas dañadas en el cable (mismo criterio que NMEA)."""
    c = 0
    for ch in texto:
        c ^= ord(ch)
    return c


class Ingestor:
    def __init__(self, con, inicio_sim):
        self.con = con
        self.inicio_sim = inicio_sim
        self.eventos = deque(maxlen=400)      # últimos eventos para la pantalla (no es la base)
        self.seq = 0
        self.lineas_ok = self.lineas_descartadas = self.ciclos_incompletos = 0
        self.ciclos = 0
        self.ciclo = None                     # ciclo en armado: cruce + medición + corrección
        self.cruces = 0
        self.ms_primer_cruce = self.ms_ultimo_cruce = None
        self.hueco_inicio = None              # (ms, azimut al perder) mientras dura un hueco
        self.huecos = 0
        self.estado = "SIGUIENDO"
        self.azimut = 0.0
        params = dict(largo_m=LARGO, latitud_grados=LATITUD_GRADOS, amplitud_m=AMPLITUD,
                      rango_sensores_grados=RANGO_SENSORES, umbral_suma_v=UMBRAL_SUMA,
                      ciclos_para_perdida=CICLOS_PARA_PERDIDA, pasos_barrido=PASOS_BARRIDO)
        cur = con.execute(
            "INSERT INTO sesiones (inicio_real,inicio_sim,k_h,k_m,grados_por_paso,d_s_m,parametros) VALUES (?,?,?,?,?,?,?)",
            (datetime.now().isoformat(timespec="seconds"), inicio_sim.isoformat(sep=" "),
             K_H, K_M, GRADOS_POR_PASO, D_S, json.dumps(params)))
        self.sesion_id = cur.lastrowid
        self.tramo_numero = 0
        self._abrir_tramo(0)                  # lanzamiento: primer tramo (la torre arranca alineada)
        con.commit()

    # ------------------------------------------------------------------ utilidades
    def marca(self, ms):
        return (self.inicio_sim + timedelta(milliseconds=ms)).isoformat(sep=" ", timespec="milliseconds")

    def _evento(self, ms, tipo, detalle, linea):
        self.seq += 1
        self.eventos.append(dict(seq=self.seq, timestamp=self.marca(ms) if ms is not None else "-",
                                 tipo=tipo, detalle=detalle, linea=linea))

    def _abrir_tramo(self, ms):
        self.tramo_numero += 1
        cur = self.con.execute("INSERT INTO tramos (sesion_id,numero,inicio_ms) VALUES (?,?,?)",
                               (self.sesion_id, self.tramo_numero, ms))
        self.tramo_id = cur.lastrowid
        self.tramo_inicio_ms = ms
        self.n = 0                            # regresión incremental del tramo (Welford)
        self.mx = self.my = self.sxx = self.sxy = self.syy = 0.0
        self.ms_ultimo_ciclo = ms

    def horas_tramo(self):
        return (self.ms_ultimo_ciclo - self.tramo_inicio_ms) / 3.6e6

    def omega(self, min_horas=0.0):
        """Ω del tramo actual: pendiente de azimut [°] contra tiempo [h], con su error estándar.
        min_horas: duración mínima del tramo para informar un valor (con pocos minutos no tiene sentido)."""
        if self.n < 10 or self.sxx <= 0 or self.horas_tramo() < min_horas:
            return None, None
        pend = self.sxy / self.sxx
        sse = max(self.syy - pend * self.sxy, 0.0)
        return pend, math.sqrt(sse / (self.n - 2) / self.sxx)

    # ------------------------------------------------------------------ entrada: una línea del cable
    def procesar(self, linea):
        # 1) Verificar integridad: "CUERPO*HH", donde HH es el checksum en hexadecimal
        try:
            cuerpo, chk = linea.rsplit("*", 1)
            if int(chk, 16) != checksum(cuerpo):
                raise ValueError("checksum")
            campos = cuerpo.split(";")
            tipo = campos[0]
            if CAMPOS_POR_TIPO.get(tipo) != len(campos):
                raise ValueError("formato")
            ms = int(campos[1])
            valores = [float(v) for v in campos[2:]]
        except (ValueError, IndexError):
            self.lineas_descartadas += 1
            self._evento(None, "Descartada", "Línea dañada: no pasa la verificación", linea)
            return
        self.lineas_ok += 1
        getattr(self, "_" + tipo)(ms, valores, linea)   # despacha según el tipo de línea

    # ------------------------------------------------------------------ un método por tipo de línea
    def _C(self, ms, v, linea):
        self.cruces += 1
        if self.ms_primer_cruce is None:
            self.ms_primer_cruce = ms
        self.ms_ultimo_cruce = ms
        hacia_torre = v[0] > 0
        if hacia_torre:                       # empieza un ciclo de medición
            if self.ciclo is not None:
                self.ciclos_incompletos += 1  # el ciclo anterior quedó sin cerrar (se perdió una línea)
            self.ciclo = {"ms": ms}
        self._evento(ms, "Cruce", "hacia la torre" if hacia_torre else "hacia el lado opuesto", linea)

    def _M(self, ms, v, linea):
        if self.ciclo is None:
            self.ciclos_incompletos += 1
            return
        self.ciclo["m"] = v                   # línea base, pico ida, pico regreso, V_OE, V_d, suma
        self._evento(ms, "Medición", f"V_OE = {v[3]:+.4f} V · suma = {v[5]:.2f} V", linea)

    def _K(self, ms, v, linea):
        if self.ciclo is None or "m" not in self.ciclo:
            self.ciclos_incompletos += 1
            self.ciclo = None
            return
        incremento, pasos, pasos_acum, acum = v[0], int(v[1]), int(v[2]), v[3]
        # azimut del tramo = pasos ejecutados + fracción pendiente, en grados (4.4)
        self.azimut = (pasos_acum + acum) * GRADOS_POR_PASO
        m = self.ciclo["m"]
        self.con.execute(
            "INSERT INTO datos_crudos (tramo_id,ms_cruce,timestamp_sim,linea_base,pico_ida,pico_regreso,v_oe,v_d,suma,"
            "incremento,pasos,pasos_acumulados,acumulador,azimut_deg) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (self.tramo_id, self.ciclo["ms"], self.marca(self.ciclo["ms"]), *m, incremento, pasos, pasos_acum, acum, self.azimut))
        self.ciclos += 1
        # regresión incremental: x = horas desde el inicio del tramo, y = azimut [°]
        xh = (self.ciclo["ms"] - self.tramo_inicio_ms) / 3.6e6
        self.n += 1
        dx = xh - self.mx
        self.mx += dx / self.n
        dy = self.azimut - self.my
        self.my += dy / self.n
        self.sxx += dx * (xh - self.mx)
        self.sxy += dx * (self.azimut - self.my)
        self.syy += dy * (self.azimut - self.my)
        self.ms_ultimo_ciclo = self.ciclo["ms"]
        self.ciclo = None
        detalle = (f"{pasos:+d} paso{'s' if abs(pasos) != 1 else ''}" if pasos else "sin movimiento") + f" · azimut {self.azimut:.2f}°"
        self._evento(ms, "Corrección", detalle, linea)

    def _P(self, ms, v, linea):
        self.hueco_inicio = (ms, (v[0] + v[1]) * GRADOS_POR_PASO)
        self.estado = "BARRIDO"
        self.ciclo = None
        self._evento(ms, "Pérdida de referencia", "los sensores no ven el imán: empieza el barrido", linea)

    def _R(self, ms, v, linea):
        """Cierre del hueco: tres escrituras en UNA transacción (4.6, 'Relaciones entre tablas')."""
        inicio_ms, az_perdida = self.hueco_inicio if self.hueco_inicio else (None, None)
        omega, se = self.omega()
        tramo_anterior = self.tramo_id
        self.con.commit()                     # cierra lo pendiente para que el bloque sea independiente
        with self.con:                        # todo o nada: si se corta la luz, no queda a medias
            self.con.execute("UPDATE tramos SET fin_ms=?, azimut_final_deg=?, omega_deg_h=?, omega_se_deg_h=?, n_ciclos=? WHERE id=?",
                             (inicio_ms if inicio_ms is not None else ms, az_perdida, omega, se, self.n, tramo_anterior))
            self._abrir_tramo(ms)
            self.con.execute("INSERT INTO huecos (tramo_cierra_id,tramo_abre_id,inicio_ms,fin_ms,pasos_barrido,azimut_al_perder_deg)"
                             " VALUES (?,?,?,?,?,?)", (tramo_anterior, self.tramo_id, inicio_ms, ms, int(v[0]), az_perdida))
        self.huecos += 1
        self.hueco_inicio = None
        self.estado = "SIGUIENDO"
        self.azimut = 0.0
        self._evento(ms, "Recuperación", f"barrido de {int(v[0])} pasos ({v[0] * GRADOS_POR_PASO:.1f}°) · tramo {self.tramo_numero}", linea)


# ==============================================================================================
# SIMULADOR: une las tres partes y avanza el tiempo
# ==============================================================================================
class Simulador:
    def __init__(self, ruta_db):
        self.ruta_db = ruta_db
        self.lock = threading.Lock()          # protege el estado entre este hilo y el servidor web
        self.velocidad = 1
        self.pausado = False
        self.cable_ruidoso = False
        self.reiniciar_pendiente = False
        self.perturbar_pendiente = False
        self.con = None
        self.ingestor = None
        self.arduino = None
        self.t = 0.0

    def _abrir_db(self):
        # Versión del esquema guardada en la propia base (PRAGMA user_version). Si el archivo es de una
        # versión anterior, no se mezcla: se renombra como respaldo y se crea una base nueva.
        if os.path.exists(self.ruta_db):
            con = sqlite3.connect(self.ruta_db)
            version = con.execute("PRAGMA user_version").fetchone()[0]
            con.close()
            if version != VERSION_ESQUEMA:
                respaldo = self.ruta_db.replace(".db", f"_respaldo_v{version}_{datetime.now():%Y%m%d_%H%M%S}.db")
                for sufijo in ("", "-wal", "-shm"):
                    if os.path.exists(self.ruta_db + sufijo):
                        os.replace(self.ruta_db + sufijo, respaldo + sufijo)
                print(f"Base de una versión anterior (esquema {version}): guardada como {os.path.basename(respaldo)}")
        self.con = sqlite3.connect(self.ruta_db)
        self.con.execute("PRAGMA journal_mode=WAL")      # lectores no bloquean al escritor (4.6)
        self.con.execute("PRAGMA foreign_keys=ON")       # SQLite lo trae apagado y es por conexión (4.6)
        self.con.executescript("""
            CREATE TABLE IF NOT EXISTS sesiones (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                inicio_real TEXT NOT NULL, inicio_sim TEXT NOT NULL,
                k_h REAL NOT NULL, k_m REAL NOT NULL, grados_por_paso REAL NOT NULL, d_s_m REAL NOT NULL,
                parametros TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS tramos (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                sesion_id INTEGER NOT NULL REFERENCES sesiones(id),
                numero INTEGER NOT NULL,
                inicio_ms INTEGER NOT NULL, fin_ms INTEGER,          -- fin vacío = tramo abierto
                azimut_final_deg REAL, omega_deg_h REAL, omega_se_deg_h REAL, n_ciclos INTEGER,
                CHECK (fin_ms IS NULL OR fin_ms >= inicio_ms));
            CREATE TABLE IF NOT EXISTS datos_crudos (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                tramo_id INTEGER NOT NULL REFERENCES tramos(id),
                ms_cruce INTEGER NOT NULL, timestamp_sim TEXT NOT NULL,
                linea_base REAL, pico_ida REAL, pico_regreso REAL, v_oe REAL, v_d REAL, suma REAL,
                incremento REAL, pasos INTEGER, pasos_acumulados INTEGER, acumulador REAL, azimut_deg REAL);
            CREATE TABLE IF NOT EXISTS huecos (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                tramo_cierra_id INTEGER NOT NULL REFERENCES tramos(id),
                tramo_abre_id   INTEGER NOT NULL REFERENCES tramos(id),
                inicio_ms INTEGER, fin_ms INTEGER NOT NULL,
                pasos_barrido INTEGER, azimut_al_perder_deg REAL);
            CREATE INDEX IF NOT EXISTS ix_crudos_tramo ON datos_crudos(tramo_id);
        """)
        self.con.execute(f"PRAGMA user_version = {VERSION_ESQUEMA}")   # marca la versión del esquema
        self.con.commit()

    def _nueva_sesion(self):
        self.t, self.resto = 0.0, 0.0
        self.x, self.y, self.vx, self.vy = AMPLITUD, 0.0, 0.0, 0.0   # se suelta desde x = AMPLITUD
        self.az_plano_real = 0.0
        self.inicio_sim = datetime.now().replace(microsecond=0)
        self.agenda = []                     # tareas del Arduino pendientes: (instante, función)
        self.ingestor = Ingestor(self.con, self.inicio_sim)
        self.arduino = ArduinoVirtual(self._cable, self._programar)

    def _programar(self, t, funcion):
        self.agenda.append((t, funcion))     # los instantes llegan en orden creciente

    def _cable(self, cuerpo):
        """El 'cable serie': agrega el checksum y, si está activado el ruido, a veces daña la línea.
        En el sistema real, en este punto el ingestor leería del puerto USB con pyserial."""
        linea = f"{cuerpo}*{checksum(cuerpo):02X}"
        if self.cable_ruidoso and random.random() < PROB_CORRUPCION:
            i = random.randrange(len(linea) - 3)
            linea = linea[:i] + random.choice("0123456789;#") + linea[i + 1:]   # un carácter cambiado
        self.ingestor.procesar(linea)

    def _avanzar(self, dt_sim):
        self.resto += dt_sim
        pasos = int(self.resto / DT)
        self.resto -= pasos * DT
        x, y, vx, vy, t = self.x, self.y, self.vx, self.vy, self.t
        for _ in range(pasos):
            rv_ant = x * vx + y * vy               # r·v: + se aleja del centro, - se acerca
            x0, y0, vx0, vy0 = x, y, vx, vy
            x, y, vx, vy = paso_rk4(x, y, vx, vy)
            t += DT
            rv = x * vx + y * vy
            if rv_ant < 0.0 <= rv:                  # mínimo de distancia al centro = cruce
                f = rv_ant / (rv_ant - rv)          # interpolación dentro del paso
                vxc, vyc = vx0 + f * (vx - vx0), vy0 + f * (vy - vy0)
                rumbo = math.degrees(math.atan2(vyc, vxc)) % 360.0       # hacia dónde va la esfera
                self.az_plano_real = rumbo % 180.0                       # el plano (verdad del mundo)
                self.arduino.al_detectar_cruce(t - DT + f * DT, self.az_plano_real, rumbo)
            while self.agenda and self.agenda[0][0] <= t:   # ejecuta lo que el Arduino tenía agendado
                self.agenda.pop(0)[1]()
        self.x, self.y, self.vx, self.vy, self.t = x, y, vx, vy, t
        self.con.commit()                           # todas las filas del intervalo en un solo bloque

    def _perturbar(self):
        """Gira el plano de golpe (como un empujón): fuerza la pérdida de referencia."""
        a = math.radians(PERTURBACION_GRADOS)
        c, s = math.cos(a), math.sin(a)
        self.x, self.y = c * self.x - s * self.y, s * self.x + c * self.y
        self.vx, self.vy = c * self.vx - s * self.vy, s * self.vx + c * self.vy
        self.ingestor._evento(round(self.t * 1000), "Simulación", f"perturbación: el plano gira {PERTURBACION_GRADOS:.0f}° de golpe", "")

    def run(self):
        self._abrir_db()
        with self.lock:
            self._nueva_sesion()
        t_real_ant = time.perf_counter()
        while True:
            time.sleep(1 / 60)
            ahora = time.perf_counter()
            dt_real = min(ahora - t_real_ant, 0.1)
            t_real_ant = ahora
            with self.lock:
                if self.reiniciar_pendiente:
                    self.reiniciar_pendiente = False
                    self._nueva_sesion()
                if self.perturbar_pendiente:
                    self.perturbar_pendiente = False
                    self._perturbar()
                if not self.pausado:
                    self._avanzar(dt_real * self.velocidad)    # tiempo simulado = real × factor

    def estado(self):
        with self.lock:
            ing, ard = self.ingestor, self.arduino
            if ing is None:
                return dict(listo=False)
            omega, se = ing.omega(MIN_HORAS_OMEGA)
            faltan_min = max(0.0, (MIN_HORAS_OMEGA - ing.horas_tramo()) * 60)
            periodo = None
            if ing.cruces >= 3:
                periodo = 2.0 * (ing.ms_ultimo_cruce - ing.ms_primer_cruce) / 1000.0 / (ing.cruces - 1)
            return dict(
                listo=True, sesion_id=ing.sesion_id, t_sim=self.t,
                timestamp_sim=(self.inicio_sim + timedelta(seconds=self.t)).isoformat(sep=" ", timespec="milliseconds"),
                x=self.x, y=self.y, vx=self.vx, vy=self.vy, largo=LARGO, d_s=D_S,
                velocidad=self.velocidad, pausado=self.pausado, cable_ruidoso=self.cable_ruidoso,
                az_torre=ard.az_torre, az_plano_real=self.az_plano_real,
                desvio_real=desvio(self.az_plano_real, ard.az_torre),
                estado=ing.estado, tramo=ing.tramo_numero, ciclos=ing.ciclos, azimut=ing.azimut,
                periodo=periodo, omega=omega, omega_se=se, omega_teorica=OMEGA_TEORICA, omega_faltan_min=faltan_min,
                cumplimiento=None if omega is None else 100.0 * (1.0 - abs(omega - OMEGA_TEORICA) / OMEGA_TEORICA),
                huecos=ing.huecos, lineas_ok=ing.lineas_ok, lineas_descartadas=ing.lineas_descartadas,
                ciclos_incompletos=ing.ciclos_incompletos, resolucion=GRADOS_POR_PASO)

    def eventos(self, desde):
        with self.lock:
            if self.ingestor is None:
                return []
            nuevos = [e for e in self.ingestor.eventos if e["seq"] > desde]
            return nuevos[-60:]                       # como máximo los 60 más recientes


# ==============================================================================================
# SERVIDOR WEB
# ==============================================================================================
app = Flask(__name__, static_folder=os.path.join(BASE, "static"), static_url_path="/static")
sim = Simulador(RUTA_DB)


@app.route("/")
def inicio():
    return send_from_directory(BASE, "index.html")


@app.route("/api/estado")
def api_estado():
    return jsonify(sim.estado())


@app.route("/api/eventos")
def api_eventos():
    return jsonify(sim.eventos(request.args.get("desde", default=0, type=int)))


@app.route("/api/control", methods=["POST"])
def api_control():
    datos = request.get_json(force=True, silent=True) or {}
    with sim.lock:
        if datos.get("velocidad") in VELOCIDADES:
            sim.velocidad = datos["velocidad"]
        if "pausado" in datos:
            sim.pausado = bool(datos["pausado"])
        if "cable_ruidoso" in datos:
            sim.cable_ruidoso = bool(datos["cable_ruidoso"])
        if datos.get("reiniciar"):
            sim.reiniciar_pendiente = True
        if datos.get("perturbar"):
            sim.perturbar_pendiente = True
    return jsonify(ok=True)


if __name__ == "__main__":
    threading.Thread(target=sim.run, daemon=True).start()
    print(f"Precesión teórica a {LATITUD_GRADOS}°: {OMEGA_TEORICA:.3f} °/h | período: {2 * math.pi / math.sqrt(W02):.3f} s")
    print("Abrí http://localhost:5000 en el navegador")
    app.run(host="127.0.0.1", port=5000, threaded=True)
