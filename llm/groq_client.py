"""Groq API client (free tier — no credit card required).

Free tier limits (as of 2025):
  - 30 requests/minute
  - 500,000 tokens/day (per model)
  - No billing required

Recommended free models:
  - llama-3.3-70b-versatile   (高品質・推奨)
  - llama-3.1-8b-instant      (高速・軽量)
  - mixtral-8x7b-32768        (長文コンテキスト)

Setup:
  pip install groq
  Get free API key: https://console.groq.com (Googleアカウントで登録可)
"""
import os
from typing import Optional

try:
    from groq import Groq
    _GROQ_AVAILABLE = True
except ImportError:
    _GROQ_AVAILABLE = False


class GroqClient:
    def __init__(self, api_key: Optional[str] = None, model: str = "llama-3.3-70b-versatile"):
        if not _GROQ_AVAILABLE:
            raise ImportError("Run: pip install groq")

        key = api_key or os.environ.get("GROQ_API_KEY")
        if not key:
            raise ValueError(
                "Groq API key required.\n"
                "1. Get free key: https://console.groq.com\n"
                "2. Set: os.environ['GROQ_API_KEY'] = 'your-key'"
            )
        self.client = Groq(api_key=key)
        self.model = model

    def generate(self, prompt: str, temperature: float = 0.3, max_tokens: int = 4096) -> str:
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[{"role": "user", "content": prompt}],
            temperature=temperature,
            max_tokens=max_tokens,
        )
        return response.choices[0].message.content
