FROM debian:12

RUN apt-get update -y && apt-get install -y cmake git g++ libtbb-dev libmpfr-dev libgmp-dev ninja-build python3-dev

ARG UID=1000
ARG GID=$UID

RUN apt-get install -y sudo && \
    addgroup --gid $GID nonroot && \
    adduser --uid $UID --gid $GID --disabled-password --gecos "" nonroot && \
    echo 'nonroot ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers

USER nonroot