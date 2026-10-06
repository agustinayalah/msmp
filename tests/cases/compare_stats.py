#!/usr/bin/env python3
"""
Compara estadísticamente la salida de sample_stats de dos corridas (ej. ms
vs mspar) para decidir si sus distribuciones son equivalentes.

Uso:
    python3 compare_stats.py case1_ms_stats.txt case1_mspar_stats.txt

Entrada: dos archivos con el formato que imprime sample_stats, una línea
por réplica:
    pi:  <val>  ss: <val>  D: <val>  thetaH: <val>  H: <val>

Para cada uno de los 5 estadísticos se compara la distribución empírica de
ambos archivos con Kolmogorov-Smirnov de dos muestras (test principal, no
asume normalidad) y Mann-Whitney U (sensible a diferencias de mediana). Los
5 p-valores se corrigen por comparaciones múltiples (Benjamini-Hochberg).

No compara réplica por réplica (no tiene sentido si cada programa usa
streams de números aleatorios distintos): compara si vienen de la misma
distribución.
"""
import re
import sys

import numpy as np
from scipy import stats

LINE_RE = re.compile(
    r"pi:\t([\-0-9.eE]+)\tss:\t(\d+)\tD:\t([\-0-9.eEnNaA]+)\tthetaH:\t([\-0-9.eE]+)\tH:\t([\-0-9.eE]+)"
)
STATS = ["pi", "ss", "D", "thetaH", "H"]


def parse(path):
    cols = {k: [] for k in STATS}
    with open(path) as f:
        for line in f:
            m = LINE_RE.match(line)
            if not m:
                continue
            pi, ss, D, thetaH, H = m.groups()
            cols["pi"].append(float(pi))
            cols["ss"].append(float(ss))
            cols["D"].append(float(D) if D.lower() != "nan" else np.nan)
            cols["thetaH"].append(float(thetaH))
            cols["H"].append(float(H))
    return {k: np.array(v)[~np.isnan(v)] if (v := cols[k]) else np.array([]) for k in cols}


def bh_adjust(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit(f"uso: {sys.argv[0]} <stats_programaA.txt> <stats_programaB.txt>")

    path_a, path_b = sys.argv[1], sys.argv[2]
    a = parse(path_a)
    b = parse(path_b)

    results = []
    for stat_name in STATS:
        xa, xb = a[stat_name], b[stat_name]
        ks_d, ks_p = stats.ks_2samp(xa, xb)
        mw_u, mw_p = stats.mannwhitneyu(xa, xb, alternative="two-sided")
        results.append({
            "stat": stat_name, "n_a": len(xa), "n_b": len(xb),
            "mean_a": xa.mean(), "mean_b": xb.mean(),
            "ks_d": ks_d, "ks_p": ks_p, "mw_p": mw_p,
        })

    ks_p_adj = bh_adjust([r["ks_p"] for r in results])
    for r, p_adj in zip(results, ks_p_adj):
        r["ks_p_adj"] = p_adj

    print(f"A = {path_a}   B = {path_b}\n")
    header = f"{'estad.':8} {'n_A':>8} {'n_B':>8} {'media_A':>14} {'media_B':>14} {'KS_D':>8} {'KS_p_adj':>10}  veredicto"
    print(header)
    print("-" * len(header))

    all_equivalent = True
    for r in results:
        equivalent = r["ks_p_adj"] > 0.05
        all_equivalent &= equivalent
        veredicto = "equivalen" if equivalent else "DIFIEREN"
        print(f"{r['stat']:8} {r['n_a']:>8} {r['n_b']:>8} {r['mean_a']:>14.4f} {r['mean_b']:>14.4f} "
              f"{r['ks_d']:>8.4f} {r['ks_p_adj']:>10.4g}  {veredicto}")

    print()
    if all_equivalent:
        print("=> No se detectan diferencias estadísticamente significativas (p ajustado > 0.05 "
              "en los 5 estadísticos). Los dos programas son consistentes entre sí.")
    else:
        peores = [r["stat"] for r in results if r["ks_p_adj"] <= 0.05]
        print(f"=> Hay diferencia estadísticamente significativa en: {', '.join(peores)}. "
              "Revisar esos estadísticos antes de concluir equivalencia.")


if __name__ == "__main__":
    main()
