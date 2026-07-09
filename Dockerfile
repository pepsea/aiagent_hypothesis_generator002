FROM python:3.11-slim

WORKDIR /app

# matplotlib(Agg)のビルド・実行に必要な最小限のシステムライブラリ
RUN apt-get update && apt-get install -y --no-install-recommends \
    libfreetype6 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt webapp/requirements_webapp.txt ./
RUN pip install --no-cache-dir -r requirements.txt -r requirements_webapp.txt

COPY . .

ENV HOST=0.0.0.0 \
    PORT=5000 \
    FLASK_DEBUG=0 \
    PYTHONUNBUFFERED=1

EXPOSE 5000

CMD ["python", "webapp/app.py"]
