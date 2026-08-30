# Skill: Gene–Disease Network Overlap Evaluator App

## 目的

このスキルは、`aiagent_hypothesis_generator002` の「NW重複（Network–Disease Gene Overlap）」解析を
独立したアプリとして別プロジェクトに移植するための実装ガイドです。

**インプット**: 疾患名 + 遺伝子リスト  
**解析**: 遺伝子ごとにPPIパートナーを取得し、疾患のOT上位遺伝子群との重複を加重スコアで評価  
**アウトプット**: 遺伝子ごとのスコア・重複遺伝子リストをテーブル表示

---

## 移植する解析ロジック（参照元ファイル）

| ファイル | 役割 |
|---|---|
| `collectors/opentargets.py` | 疾患上位遺伝子取得（`get_disease_top_genes`）+ 遺伝子ID解決（`SEARCH_QUERY`） |
| `webapp/app.py` | NW重複計算ロジック（`network_disease_overlap` ブロック） |

---

## 実装手順

### Step 1. プロジェクト構成

```
nw_overlap_app/
├── app.py              # Flask アプリ本体
├── collectors/
│   ├── __init__.py
│   ├── opentargets.py  # 元ファイルからそのままコピー
│   └── ppi.py          # 元プロジェクトの PPI 収集ロジック
├── templates/
│   └── index.html
└── requirements.txt
```

`requirements.txt`:
```
flask
requests
networkx
```

---

### Step 2. 疾患上位遺伝子取得

`collectors/opentargets.py` から以下をそのままコピーして使う：

- `get_disease_top_genes(disease_id, top_n=100)` — OTスコア上位遺伝子を返す
- `SEARCH_QUERY` + 検索ロジック — 疾患名 → EFO ID 変換に使う

疾患名 → EFO ID 変換の関数を追加：

```python
def resolve_disease_id(disease_name: str) -> tuple[str, str]:
    """疾患名 → (efo_id, label) を返す。見つからなければ (None, disease_name)。"""
    r = requests.post(OT_API, json={
        "query": SEARCH_QUERY,
        "variables": {"q": disease_name, "entity": ["disease"]}
    }, timeout=20)
    r.raise_for_status()
    hits = [h for h in r.json().get("data", {}).get("search", {}).get("hits", [])
            if h.get("entity") == "disease"]
    if not hits:
        return None, disease_name
    exact = [h for h in hits if h.get("name", "").lower() == disease_name.lower()]
    best = exact[0] if exact else hits[0]
    return best["id"], best["name"]
```

---

### Step 3. PPI パートナー取得

元プロジェクトの `webapp/app.py` の PPI 収集ブロックを関数化する。
SIGNOR / STRING / BioGRID を使う最小版：

```python
# collectors/ppi.py
import requests, networkx as nx

SIGNOR_URL = "https://signor.uniroma2.it/getData.php"
STRING_URL = "https://string-db.org/api/json/network"
BIOGRID_URL = "https://webservice.thebiogrid.org/interactions/"

def get_ppi_partners(gene: str, biogrid_key: str = "", top_n: int = 30) -> list[str]:
    """SIGNOR + STRING + BioGRID から PPI パートナーを収集して上位 top_n を返す。"""
    partners = set()

    # SIGNOR
    try:
        r = requests.get(SIGNOR_URL, params={"organism": "9606", "format": "json"}, timeout=15)
        for row in r.json():
            a, b = row.get("ENTITYA", ""), row.get("ENTITYB", "")
            if a.upper() == gene.upper(): partners.add(b)
            elif b.upper() == gene.upper(): partners.add(a)
    except Exception:
        pass

    # STRING (requires Ensembl ID or gene name)
    try:
        r = requests.get(STRING_URL, params={
            "identifiers": gene, "species": 9606, "limit": top_n, "caller_identity": "nw_overlap_app"
        }, timeout=15)
        for row in r.json():
            a, b = row.get("preferredName_A", ""), row.get("preferredName_B", "")
            if a.upper() == gene.upper(): partners.add(b)
            elif b.upper() == gene.upper(): partners.add(a)
    except Exception:
        pass

    # BioGRID (API key 必要)
    if biogrid_key:
        try:
            r = requests.get(BIOGRID_URL, params={
                "geneList": gene, "taxId": 9606, "accesskey": biogrid_key,
                "format": "json", "max": top_n, "interSpeciesExcluded": "true"
            }, timeout=15)
            for v in r.json().values():
                a = v.get("OFFICIAL_SYMBOL_A", "")
                b = v.get("OFFICIAL_SYMBOL_B", "")
                if a.upper() == gene.upper(): partners.add(b)
                elif b.upper() == gene.upper(): partners.add(a)
        except Exception:
            pass

    partners.discard(gene.upper())
    return list(partners)[:top_n]
```

---

### Step 4. NW重複スコア計算

元プロジェクトの計算ロジックをそのまま関数化する：

```python
def calc_network_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_genes: list[dict],   # [{"symbol": str, "score": float}, ...]
) -> dict:
    """
    NW重複スコアを計算する。

    weighted_score = Σ(OTスコア of overlapping genes + target自身のスコア) / Σ(全疾患遺伝子OTスコア)

    Returns:
        weighted_score, simple_ratio, overlap_count, overlapping_genes,
        target_self (target が疾患遺伝子リストにあれば symbol), target_self_score
    """
    _total_ot_score = sum(g.get("score", 0) for g in disease_genes)
    _ppi_set = {p.upper() for p in ppi_partners}
    _dg = disease_genes  # [{"symbol", "score"}, ...]

    # ターゲット自身が疾患遺伝子リストにあるか確認
    _self_entry = next((g for g in _dg if g.get("symbol", "").upper() == gene.upper()), None)
    _self_score = _self_entry.get("score") if _self_entry else None

    # PPIパートナーと疾患遺伝子の重複（自身は除く）
    _overlap = [
        g for g in _dg
        if g.get("symbol", "").upper() in _ppi_set
        and g.get("symbol", "").upper() != gene.upper()
    ]

    _weighted = sum(g.get("score", 0) for g in _overlap)
    if _self_score is not None:
        _weighted += _self_score

    _simple_n = len(_overlap) + (1 if _self_entry else 0)
    _simple_ratio = round(_simple_n / max(1, len(_dg)), 3)
    _weighted_score = round(_weighted / _total_ot_score, 3) if _total_ot_score > 0 else 0.0

    return {
        "weighted_score":   _weighted_score,
        "simple_ratio":     _simple_ratio,
        "overlap_count":    len(_overlap),
        "disease_gene_count": len(_dg),
        "ppi_partner_count":  len(ppi_partners),
        "target_self":        gene if _self_entry else None,
        "target_self_score":  _self_score,
        "overlapping_genes":  sorted(_overlap, key=lambda g: g.get("score", 0), reverse=True),
    }
```

---

### Step 5. Flask アプリ本体（app.py）

```python
from flask import Flask, request, render_template, jsonify
from collectors.opentargets import resolve_disease_id, get_disease_top_genes
from collectors.ppi import get_ppi_partners
from nw_overlap import calc_network_overlap
import concurrent.futures

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/api/analyze", methods=["POST"])
def analyze():
    data = request.get_json()
    disease_name = (data.get("disease") or "").strip()
    genes = [g.strip() for g in (data.get("genes") or []) if g.strip()]

    if not disease_name or not genes:
        return jsonify({"error": "disease and genes are required"}), 400

    # 疾患ID解決
    disease_id, disease_label = resolve_disease_id(disease_name)
    if not disease_id:
        return jsonify({"error": f"Disease not found: {disease_name}"}), 404

    # 疾患上位遺伝子取得（上位100件）
    disease_genes = get_disease_top_genes(disease_id, top_n=100)

    # 各遺伝子を並列処理
    results = []
    biogrid_key = app.config.get("BIOGRID_KEY", "")

    def process_gene(gene):
        partners = get_ppi_partners(gene, biogrid_key=biogrid_key)
        overlap = calc_network_overlap(gene, partners, disease_genes)
        return {"gene": gene, **overlap}

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
        futures = {ex.submit(process_gene, g): g for g in genes}
        for f in concurrent.futures.as_completed(futures):
            try:
                results.append(f.result())
            except Exception as e:
                results.append({"gene": futures[f], "error": str(e)})

    results.sort(key=lambda r: r.get("weighted_score", 0), reverse=True)

    return jsonify({
        "disease_id":    disease_id,
        "disease_label": disease_label,
        "disease_gene_count": len(disease_genes),
        "results": results,
    })

if __name__ == "__main__":
    app.run(debug=True, port=5010)
```

---

### Step 6. フロントエンド（templates/index.html）

#### 入力フォーム

```html
<form id="form">
  <label>疾患名</label>
  <input id="disease" type="text" placeholder="例: Alzheimer disease" required>

  <label>遺伝子リスト（1行1遺伝子 または カンマ区切り）</label>
  <textarea id="genes" rows="6" placeholder="APP&#10;PSEN1&#10;APOE"></textarea>

  <button type="submit">解析実行</button>
</form>
<div id="result"></div>
```

#### 結果テーブル描画

```javascript
async function runAnalysis(e) {
  e.preventDefault();
  const disease = document.getElementById("disease").value.trim();
  const genesRaw = document.getElementById("genes").value;
  // 1行1遺伝子 or カンマ区切り両対応
  const genes = genesRaw.split(/[\n,]+/).map(g => g.trim()).filter(Boolean);

  document.getElementById("result").innerHTML = "<p>解析中...</p>";

  const res = await fetch("/api/analyze", {
    method: "POST",
    headers: {"Content-Type": "application/json"},
    body: JSON.stringify({disease, genes}),
  });
  const data = await res.json();
  if (data.error) { document.getElementById("result").innerHTML = `<p style="color:red">${data.error}</p>`; return; }

  let html = `<p><strong>${data.disease_label}</strong>（${data.disease_id}）— 疾患上位遺伝子 ${data.disease_gene_count} 件</p>`;
  html += `<table><thead><tr>
    <th>遺伝子</th><th>加重スコア</th><th>単純重複率</th>
    <th>重複数</th><th>ターゲット自身</th><th>重複遺伝子（上位5件）</th>
  </tr></thead><tbody>`;

  for (const r of data.results) {
    if (r.error) { html += `<tr><td>${r.gene}</td><td colspan="5" style="color:red">${r.error}</td></tr>`; continue; }
    const selfMark = r.target_self ? `★ ${r.target_self} (${r.target_self_score?.toFixed(3)})` : "—";
    const topGenes = (r.overlapping_genes || []).slice(0, 5).map(g => `${g.symbol}(${g.score.toFixed(3)})`).join(", ");
    html += `<tr>
      <td><strong>${r.gene}</strong></td>
      <td>${r.weighted_score.toFixed(3)}</td>
      <td>${(r.simple_ratio * 100).toFixed(1)}%</td>
      <td>${r.overlap_count} / ${r.ppi_partner_count} PPI</td>
      <td>${selfMark}</td>
      <td style="font-size:12px">${topGenes || "—"}</td>
    </tr>`;
  }
  html += "</tbody></table>";
  document.getElementById("result").innerHTML = html;
}
document.getElementById("form").addEventListener("submit", runAnalysis);
```

---

## 実装時の注意点

### 入力バリデーション
- 疾患名：空文字チェック → OT検索でヒットしなければ 404 返却
- 遺伝子リスト：1行1遺伝子 + カンマ区切りの両方をパース、空行・重複は除去
- 最大遺伝子数：並列処理の負荷を考慮して 20〜30 件を上限推奨

### OpenTargets `top_n` の選択
- 疾患上位遺伝子は `top_n=100` を推奨（多すぎると重複率が高止まりし識別力が落ちる）
- ニッチな疾患では遺伝子数が少なく `top_n=50` でも十分

### PPIパートナー数
- SIGNOR/STRING 合計で 20〜50 件が一般的な範囲
- BioGRID は API キーが必要。なければ SIGNOR + STRING のみでも動作する

### スコアの解釈
| 加重スコア | 解釈 |
|---|---|
| ≥ 0.3 | ネットワーク的に疾患と強く関連 |
| 0.1〜0.3 | 中程度の関連 |
| < 0.1 | 関連が薄い（または PPI 情報不足） |

---

## 元プロジェクトとの差分

| 項目 | 元プロジェクト | この独立アプリ |
|---|---|---|
| 入力 | 遺伝子1件 × 疾患1件（逐次） | 遺伝子Nリスト × 疾患1件（並列） |
| LLM | 仮説生成あり | なし |
| 表示 | SSE ストリーミング | REST API + 静的テーブル |
| 保存 | HTML スナップショット | なし（必要なら追加） |
| PPI 除外ハブ | あり（degree > 閾値で除外） | 簡略化（なし） |
