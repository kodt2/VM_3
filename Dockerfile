FROM python:3.12-slim

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir fastapi==0.115.6 uvicorn[standard]==0.34.0

COPY solver.cpp main.py index.html ./
RUN g++ -std=c++17 -O2 -Wall -Wextra -pedantic solver.cpp -o solver

EXPOSE 8000
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
