"""Ollama local LLM client (completely free, runs locally).

Setup:
  1. Install: https://ollama.com
  2. Pull model: ollama pull llama3.1
  3. Start server: ollama serve  (or it starts automatically)
"""
import requests
from typing import Optional


class OllamaClient:
    def __init__(self, model: str = "llama3.1", base_url: str = "http://localhost:11434"):
        self.model = model
        self.base_url = base_url

    def generate(self, prompt: str, temperature: float = 0.3, max_tokens: int = 4096) -> str:
        r = requests.post(f"{self.base_url}/api/generate", json={
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {
                "temperature": temperature,
                "num_predict": max_tokens,
            }
        }, timeout=120)
        r.raise_for_status()
        return r.json().get("response", "")

    def is_available(self) -> bool:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=3)
            return r.status_code == 200
        except Exception:
            return False
