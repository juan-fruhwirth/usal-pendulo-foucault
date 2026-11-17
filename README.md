# Simulador del péndulo de Foucault — USAL Pilar

Simulador del sistema de control y medición de la precesión del péndulo de Foucault de la Universidad del Salvador, sede Pilar. Reproduce el péndulo con la rotación de la Tierra, un Arduino virtual con la lógica de control del diseño y el software que correría en la Raspberry Pi (ingesta, base de datos SQLite, cálculo de la precesión y tablero 3D). Proyecto Final de Ingeniería en Informática, USAL, 2026.

## Requisitos

- Python 3.10 o superior
- Flask (única dependencia)
- Un navegador web. No necesita conexión a internet.

## Instalación

```bash
git clone https://github.com/<usuario>/usal-pendulo-foucault.git
cd usal-pendulo-foucault
python -m pip install -r requirements.txt
```

En Windows, si `python` no se reconoce, usar `py` en su lugar.

## Ejecución

```bash
python simulador.py
```

Abrir <http://localhost:5000> en el navegador. Para detenerlo: Ctrl+C en la terminal.

En pantalla:

- **Velocidad** de 1x a 1000x.
- **Perturbar plano**: gira el plano 25° de golpe; los sensores pierden el imán, la torre barre hasta reencontrarlo y empieza un tramo nuevo.
- **Cable con ruido**: daña el 2 % de los mensajes; el sistema los detecta con la suma de verificación y los descarta.

La base de datos se guarda en `pendulo.db` y se puede abrir con DB Browser for SQLite.

## Corridas del capítulo 7

Los programas de `scripts/` reproducen los resultados del documento, con semillas aleatorias fijas:

```bash
python scripts/corridas_resultados.py    # referencia, perturbación y cable con ruido (ruido mínimo)
python scripts/corridas_ruido_2_7.py     # referencia y perturbación con 2,7° de ruido por lectura
python scripts/analisis_ruido.py         # tiempo de registro según el nivel de ruido
```

## Estructura

| Archivo | Contenido |
|---|---|
| `simulador.py` | Física del péndulo, Arduino virtual, ingestor, base de datos y servidor web |
| `index.html` | Tablero: péndulo en 3D, indicadores y registro de eventos |
| `static/three.min.js` | Biblioteca Three.js incluida para funcionar sin internet |
| `scripts/` | Programas de las corridas del capítulo 7 |

## Autor

Juan Andrés Fruhwirth — Universidad del Salvador, Facultad de Ingeniería.

## Licencia

MIT. Ver `LICENSE`. Three.js se distribuye con su propia licencia MIT.
