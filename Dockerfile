FROM ubuntu:22.04 as build

RUN apt-get update -y && apt-get install -y cmake git g++ libtbb-dev libmpfr-dev libgmp-dev

RUN git clone --depth 1 https://github.com/TLCFEM/vpmr.git

WORKDIR /vpmr

RUN git submodule update --init --recursive
RUN cd eigen && git apply ../patch_size.patch && cd ..
RUN sed -i 's/3.24/3.13/g' CMakeLists.txt
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

FROM ubuntu:22.04 as runtime

RUN apt-get update -y && apt-get install -y libtbb12 libmpfr6 libgmp10 && apt-get clean -y

COPY --from=build /vpmr/build/vpmr /usr/local/bin/vpmr

RUN vpmr -h

ENTRYPOINT ["vpmr"]
