# 岩芯 CT 一维持久同源分析 —— 纯后端服务镜像
FROM python:3.11-slim

# 容器内监听端口可由 PORT 覆盖；宿主机映射端口在 docker-compose.yml 中配置
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8000

WORKDIR /srv/app

# 先装依赖，利用层缓存
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# 拷贝应用代码
COPY app ./app

# 以非 root 用户运行
RUN useradd --create-home --uid 10001 appuser \
    && chown -R appuser:appuser /srv/app
USER appuser

EXPOSE 8000

# 容器级健康检查（不依赖外部工具，仅用标准库）
HEALTHCHECK --interval=15s --timeout=5s --start-period=10s --retries=3 \
    CMD python -c "import json,os,sys,urllib.request; p=os.environ.get('PORT','8000'); r=urllib.request.urlopen('http://127.0.0.1:%s/health'%p,timeout=3); sys.exit(0 if r.status==200 and json.load(r).get('status')=='ok' else 1)"

# shell 形式以便 ${PORT} 展开
CMD uvicorn app.main:app --host 0.0.0.0 --port ${PORT}
