"""Ollama local LLM client (completely free, runs locally).

Setup:
  1. Install: https://ollama.com
  2. Pull model: ollama pull llama3.1
  3. Start server: ollama serve  (or it starts automatically)
"""
import json
import requests
from typing import Optional, Callable


class OllamaClient:
    def __init__(self, model: str = "llama3.1", base_url: str = "http://localhost:11434"):
        self.model = model
        self.base_url = base_url

    def generate(
        self,
        prompt: str,
        temperature: float = 0.3,
        max_tokens: int = 4096,
        stream_callback: Optional[Callable[[str], None]] = None,
        format: Optional[str] = None,
    ) -> str:
        """Generate text. If stream_callback is provided, stream token-by-token.

        format="json" forces Ollama to emit syntactically valid JSON.
        """
        use_stream = stream_callback is not None
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": use_stream,
            "options": {
                "temperature": temperature,
                "num_predict": max_tokens,
            },
        }
        if format:
            payload["format"] = format
        r = requests.post(
            f"{self.base_url}/api/generate",
            json=payload,
            stream=use_stream,
            timeout=300,
        )
        r.raise_for_status()

        if use_stream:
            full_text = []
            for line in r.iter_lines():
                if not line:
                    continue
                chunk = json.loads(line)
                token = chunk.get("response", "")
                if token:
                    full_text.append(token)
                    stream_callback(token)
                if chunk.get("done"):
                    break
            return "".join(full_text)
        else:
            return r.json().get("response", "")

    def is_available(self) -> bool:
        try:
            r = requests.get(f"{self.base_url}/api/tags", timeout=3)
            return r.status_code == 200
        except Exception:
            return False
