"""Gemini via google-genai SDK (新SDK・AI Studio キー対応).

旧 google-generativeai は廃止予定のため新 google-genai を使用。
AI Studio の無料キーで動作確認済み。

Free tier (AI Studio キー):
  - 15 req/min, 1,000,000 tokens/day
  - https://aistudio.google.com/app/apikey

Setup:
  pip install google-genai
"""
import os
import time
import re
from typing import Optional

try:
    from google import genai
    from google.genai.errors import ClientError
    _GENAI_AVAILABLE = True
except ImportError:
    _GENAI_AVAILABLE = False


class GeminiClient:
    def __init__(self, api_key: Optional[str] = None, model: str = "gemini-2.0-flash"):
        if not _GENAI_AVAILABLE:
            raise ImportError("Run: pip install google-genai")

        key = api_key or os.environ.get("GEMINI_API_KEY")
        if not key:
            raise ValueError(
                "Gemini API key required.\n"
                "Get free key: https://aistudio.google.com/app/apikey"
            )
        self.client = genai.Client(api_key=key)
        self.model = model

    def generate(self, prompt: str, temperature: float = 0.3, max_tokens: int = 4096,
                 max_retries: int = 3) -> str:
        from google.genai import types
        config = types.GenerateContentConfig(
            temperature=temperature,
            max_output_tokens=max_tokens,
        )
        for attempt in range(max_retries):
            try:
                response = self.client.models.generate_content(
                    model=self.model,
                    contents=prompt,
                    config=config,
                )
                return response.text
            except ClientError as e:
                if "429" in str(e) or "RESOURCE_EXHAUSTED" in str(e):
                    wait = 65
                    m = re.search(r"retry in (\d+)", str(e))
                    if m:
                        wait = int(m.group(1)) + 5
                    if attempt < max_retries - 1:
                        print(f"  ⚠ レート制限。{wait}秒後にリトライ... ({attempt+1}/{max_retries})")
                        time.sleep(wait)
                    else:
                        raise RuntimeError(
                            "Gemini APIのレート制限に繰り返し達しました。\n"
                            "Groq API (llm/groq_client.py) への切り替えを推奨します。"
                        ) from e
                else:
                    raise
