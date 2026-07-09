"""Human Protein Atlas — tissue/cell/disease expression (CC BY-SA 4.0, 商用利用可).

Compact API: https://www.proteinatlas.org/{ENSG_ID}.json
    "RNA tissue specific nTPM" 等の compact JSON フィールドは、その遺伝子が
    "tissue enriched/enhanced" と判定された組織のみを返す（HPA の分類ロジック）。
    TP53 や ACTB のように広く均一に発現する遺伝子（"Low tissue specificity"）は
    この値が常に null になり、全組織のデータが欠落して見える。

    全組織・全細胞種・全がん種の完全なプロファイルは compact JSON API では
    提供されておらず、遺伝子ページ（/tissue, /single+cell, /cancer）に埋め込まれた
    棒グラフ用 JSON（$('#...').barChart([...])）にのみ含まれるため、
    これらのページを取得してスクレイピングする。
"""
import json
import re
import time

import requests
from collectors._ensembl import resolve_ensg

HPA_BASE = "https://www.proteinatlas.org"

_BAR_CHART_RE = re.compile(r"\$\('#\w+'\)\.barChart\((\[.*?\])\s*,\s*\{", re.DOTALL)
_ORGAN_RE = re.compile(r"Organ:\s*([^<]+?)(?:<br|$)")
_IHC_COUNT_RE = re.compile(
    r"High:\s*(\d+)<br>Medium:\s*(\d+)<br>Low:\s*(\d+)<br>Not detected:\s*(\d+)")


def _fetch(url: str, max_retries: int = 3) -> str | None:
    """GET url、一時的なエラーはリトライ、404 は None を返す。"""
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=20)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.text
        except Exception:
            if attempt < max_retries - 1:
                time.sleep(1 + attempt * 2)
    return None


def _first_bar_chart_after(html: str, anchor: str) -> list[dict]:
    """anchor 文字列以降で最初に現れる barChart([...]) の配列をパースして返す。"""
    idx = html.find(anchor)
    if idx == -1:
        return []
    m = _BAR_CHART_RE.search(html, idx)
    if not m:
        return []
    try:
        return json.loads(m.group(1))
    except Exception:
        return []


def _first_bar_chart(html: str) -> list[dict]:
    m = _BAR_CHART_RE.search(html)
    if not m:
        return []
    try:
        return json.loads(m.group(1))
    except Exception:
        return []


def _get_tissue_expression(ensg: str) -> list[dict]:
    """全組織の RNA 発現 (nTPM) を /tissue ページから取得する。"""
    html = _fetch(f"{HPA_BASE}/{ensg}/tissue")
    if not html:
        return []
    items = _first_bar_chart_after(html, 'href="#rna_expression"')
    out = []
    for it in items:
        tooltip = it.get("tooltip", "")
        m = _ORGAN_RE.search(tooltip)
        try:
            tpm = float(it.get("value", 0) or 0)
        except (TypeError, ValueError):
            tpm = 0.0
        out.append({
            "tissue": it.get("label", ""),
            "organ": m.group(1) if m else "",
            "level": f"{tpm:.1f} nTPM",
            "tpm": tpm,
        })
    out.sort(key=lambda x: x["tpm"], reverse=True)
    return out


def _get_single_cell_expression(ensg: str) -> list[dict]:
    """全細胞種の RNA 発現 (nCPM) を /single+cell ページから取得する。"""
    html = _fetch(f"{HPA_BASE}/{ensg}/single+cell")
    if not html:
        return []
    items = _first_bar_chart_after(html, 'href="#single_cell_type_summary"')
    out = []
    for it in items:
        try:
            ncpm = float(it.get("value", 0) or 0)
        except (TypeError, ValueError):
            ncpm = 0.0
        out.append({
            "cell_type": it.get("label", ""),
            "group": it.get("legend", ""),
            "level": f"{ncpm:.1f} nCPM",
            "ncpm": ncpm,
        })
    out.sort(key=lambda x: x["ncpm"], reverse=True)
    return out


def _get_cancer_expression(ensg: str) -> list[dict]:
    """全がん種の免疫組織化学 (IHC) 発現データを /cancer ページから取得する。"""
    html = _fetch(f"{HPA_BASE}/{ensg}/cancer")
    if not html:
        return []
    items = _first_bar_chart(html)
    out = []
    for it in items:
        tooltip = it.get("tooltip", "")
        m = _IHC_COUNT_RE.search(tooltip)
        high, medium, low, not_detected = (int(x) for x in m.groups()) if m else (0, 0, 0, 0)
        total = high + medium + low + not_detected
        out.append({
            "cancer_type": it.get("label", ""),
            "organ": it.get("legend", ""),
            "high": high, "medium": medium, "low": low, "not_detected": not_detected,
            "n_patients": total,
            "pct_high_medium": round(100 * (high + medium) / total, 1) if total else 0.0,
        })
    out.sort(key=lambda x: x["pct_high_medium"], reverse=True)
    return out


def get_expression_profile(gene_symbol: str, max_retries: int = 3) -> dict:
    """Return HPA tissue/single-cell/cancer expression and subcellular localisation.

    Returns:
        {
          "tissue_expression":     [{tissue, organ, level, tpm}],   # 全組織 RNA nTPM
          "single_cell_expression":[{cell_type, group, level, ncpm}], # 全細胞種 RNA nCPM
          "cancer_expression":     [{cancer_type, organ, high, medium, low,
                                      not_detected, n_patients, pct_high_medium}],
          "protein_tissue":        [{tissue, level}],   # protein intensity (qualitative)
          "subcellular":           [str],
          "protein_class":         [str],
          "url":                   str,
        }

    404（遺伝子が HPA に存在しない）は即座に諦めるが、タイムアウトや一時的な
    通信エラーは複数回リトライしてから諦める。
    """
    ensg = resolve_ensg(gene_symbol)
    if not ensg:
        return {"error": f"Could not resolve ENSG for {gene_symbol}"}

    # サブセルラー局在・タンパク質クラス等はコンパクト JSON API で十分
    compact_url = f"{HPA_BASE}/{ensg}.json"
    data = None
    last_err = None
    for attempt in range(max_retries):
        try:
            r = requests.get(compact_url, timeout=20)
            if r.status_code == 404:
                return {"error": f"{gene_symbol} ({ensg}) not found in HPA"}
            r.raise_for_status()
            data = r.json()
            break
        except Exception as e:
            last_err = e
            if attempt < max_retries - 1:
                time.sleep(1 + attempt * 2)
    if data is None:
        return {"error": str(last_err)}

    # Protein tissue intensity (qualitative) — compact API のみに存在
    prot_tissue = data.get("Protein tissue specific Intensity") or {}
    prot_expr = [{"tissue": t, "level": v} for t, v in prot_tissue.items()]

    # Subcellular localisation
    subcellular = []
    for k in ("Subcellular main location", "Subcellular additional location", "Subcellular location"):
        val = data.get(k)
        if isinstance(val, list):
            subcellular.extend(val)
        elif isinstance(val, str) and val:
            subcellular.append(val)
    subcellular = list(dict.fromkeys(subcellular))

    protein_class = data.get("Protein class") or []
    if isinstance(protein_class, str):
        protein_class = [protein_class]

    return {
        "tissue_expression":      _get_tissue_expression(ensg),
        "single_cell_expression": _get_single_cell_expression(ensg),
        "cancer_expression":      _get_cancer_expression(ensg),
        "protein_tissue":         prot_expr,
        "subcellular":            subcellular,
        "protein_class":          protein_class,
        "url": f"{HPA_BASE}/{ensg}",
    }
