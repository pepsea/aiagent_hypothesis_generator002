# 創薬仮説生成システム / Drug Hypothesis Generator

ローカル LLM（Ollama）と 15 以上の公開バイオデータベースを組み合わせ、
**「遺伝子 × 疾患」ごとに創薬仮説レポートを自動生成する**システムです。

- 入力: 疾患名（1件）＋ 対象遺伝子リスト（複数可）
- 処理: 公開DBからエビデンス収集 → PPIネットワーク・機能エンリッチメント解析 → LLMで仮説生成
- 出力: 遺伝子ごとの仮説レポート（PPI図・エンリッチメント・仮説本文）
- 特徴: 完全ローカル LLM（API課金なし・データ外部送信なし）。使用データは原則 CC BY 4.0（商用利用可）

3 つの使い方を用意しています。

| 形態 | ファイル | 用途 |
|---|---|---|
| **Web アプリ** | `webapp/` | ブラウザで疾患検索・遺伝子入力・リアルタイム表示（推奨） |
| ノートブック | `batch_hypothesis_generator.ipynb` | Jupyter 上でバッチ実行 |
| 単一ファイル版 | `hypothesis_generator_standalone.ipynb` | 全処理を 1 枚に収めた自己完結ノートブック |

---

## 必要環境

- Python 3.10+
- [Ollama](https://ollama.com)（ローカル LLM 実行）

```bash
# Ollama をインストール後、モデルを取得
ollama pull gemma3:4b-it-qat      # 軽量・高速（推奨）
# または
ollama pull qwen2.5:14b           # 高品質（要メモリ）
```

---

## Web アプリの起動

```bash
git clone https://github.com/pepsea/aiagent_hypothesis_generator002.git
cd aiagent_hypothesis_generator002

# 依存パッケージ
pip install -r requirements.txt
pip install -r webapp/requirements_webapp.txt

# Ollama を起動（別ターミナル）
ollama serve

# Web アプリを起動
python webapp/app.py
# → http://localhost:5000 をブラウザで開く
```

### 画面の流れ

1. **モデル選択** — Ollama にある LLM から選択
2. **疾患選択** — 疾患名を検索（OpenTargets）して確定
3. **遺伝子リスト** — HGNC シンボルを入力し検証
4. **PPI 設定** — データソース（SIGNOR / STRING / BioGRID）とクライテリア
   （STRING 信頼度閾値・エッジ最小スコア・ハブ除外閾値・最大パートナー数）を指定
5. **解析開始** — 遺伝子ごとに以下がリアルタイム表示されます
   - **取得データ** タブ: 各ソースの収集状況と内容（クリックで詳細）。
     SIGNOR/STRING/BioGRID の結果には対象遺伝子・パートナーの UniProt 機能も併記
   - **仮説** タブ: LLM 生成をストリーミング表示
   - **PPI** タブ: タンパク質相互作用ネットワーク図＋パートナー一覧＋機能リスト
   - **エンリッチメント** タブ: GO / Reactome 機能解析（ヒット遺伝子付き、ハブ遺伝子は除外）
   - **ダウンロード** タブ: Markdown レポート

### BioGRID を使う場合

BioGRID は非商用・学術利用限定ライセンスのため、API キーを環境変数で
指定したときのみ有効になります（未設定時は取得データに理由付きで
表示されます）。

```bash
# https://webservice.thebiogrid.org/ で無料登録してキーを取得
export BIOGRID_API_KEY=あなたのキー
python webapp/app.py
```

---

## データソース

PubMed / OpenTargets / UniProt / SIGNOR / STRING / BioGRID（任意） / GWAS Catalog /
ClinVar / ChEMBL / gnomAD / GTEx / Human Protein Atlas / DGIdb /
ClinicalTrials.gov / AlphaFold / g:Profiler（エンリッチメント）

---

## アーキテクチャ

```
webapp/app.py (Flask + SSE)  または  *.ipynb
        │
   pipeline.py                     ← 処理フロー制御
        │
   ┌────┼──────────┬────────────┬──────────┐
   ▼    ▼          ▼            ▼          ▼
aggregator  network      hypothesis    report    llm/ollama_client
 (収集+整形) (PPI+解析)   (プロンプト+  (Markdown  (LLM呼び出し)
                          LLM生成)      整形)
   │
   ▼
collectors/ × 18                   ← 各DBのAPIラッパー
```

詳細は [SPEC.md](SPEC.md) を参照してください。

---

## ライセンス / 注意

- 本システムが生成する仮説は **研究補助を目的とした自動生成物** であり、
  医学的・臨床的判断の根拠とするものではありません。
- 使用する公開DBはそれぞれのライセンスに従います。
