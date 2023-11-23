FROM debian as builder
RUN apt-get update && apt-get install -y cmake git g++ libtbb-dev libmpfr-dev libgmp-dev wget gtk-update-icon-cache
RUN git clone --depth 1 https://github.com/TLCFEM/vpmr.git
WORKDIR /vpmr
RUN git submodule update --init --recursive
RUN sed -i 's/3.24/3.18/g' CMakeLists.txt
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
SHELL ["/bin/bash", "-c"]
RUN echo $'version: 1 \n\
AppDir:\n\
  path: build/AppDir\n\
  app_info:\n\
    id: io.github.tlcfem.vpmr\n\
    name: vpmr\n\
    icon: vpmr\n\
    exec: usr/bin/vpmr\n\
    exec_args: $@\n\
  files:\n\
    include:\n\
    - lib64/ld-linux-x86-64.so.2\n\
    exclude:\n\
    - usr/share/man\n\
    - usr/share/doc/*/README.*\n\
    - usr/share/doc/*/changelog.*\n\
    - usr/share/doc/*/NEWS.*\n\
    - usr/share/doc/*/TODO.*\n\
  test:\n\
    fedora-30:\n\
      image: appimagecrafters/tests-env:fedora-30\n\
      command: ./AppRun -h\n\
    debian-stable:\n\
      image: appimagecrafters/tests-env:debian-stable\n\
      command: ./AppRun -h\n\
    archlinux-latest:\n\
      image: appimagecrafters/tests-env:archlinux-latest\n\
      command: ./AppRun -h\n\
    centos-7:\n\
      image: appimagecrafters/tests-env:centos-7\n\
      command: ./AppRun -h\n\
    ubuntu-xenial:\n\
      image: appimagecrafters/tests-env:ubuntu-xenial\n\
      command: ./AppRun -h\n\
AppImage:\n\
  arch: x86_64\n\
  update-information: guess\n' > AppImageBuilder.yml

RUN wget -q https://github.com/linuxdeploy/linuxdeploy/releases/download/1-alpha-20220822-1/linuxdeploy-x86_64.AppImage
RUN chmod +x linuxdeploy-x86_64.AppImage
RUN wget -q https://github.com/AppImageCrafters/appimage-builder/releases/download/v1.1.0/appimage-builder-1.1.0-x86_64.AppImage
RUN chmod +x appimage-builder-1.1.0-x86_64.AppImage

RUN ./linuxdeploy-x86_64.AppImage --appimage-extract-and-run --appdir AppDir --executable build/vpmr -i resource/vpmr.svg
RUN ./appimage-builder-1.1.0-x86_64.AppImage --appimage-extract-and-run --recipe AppImageBuilder.yml
RUN chmod +x vpmr-0.1.0-x86_64.AppImage
RUN ./vpmr-0.1.0-x86_64.AppImage --appimage-extract-and-run -h
