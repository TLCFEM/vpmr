FROM ubuntu:22.04 AS build

RUN apt-get update -y && apt-get install -y cmake git g++ libtbb-dev libmpfr-dev libgmp-dev

RUN git clone --recurse-submodules --depth 1 https://github.com/TLCFEM/vpmr.git

WORKDIR /vpmr

RUN cd eigen && git apply ../patch_size.patch && cd ..
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

FROM ubuntu:22.04 AS runtime

RUN apt-get update -y && apt-get install -y libtbb12 libmpfr6 libgmp10 && apt-get clean -y

COPY --from=build /vpmr/build/vpmr /usr/local/bin/vpmr

RUN vpmr -h

ENTRYPOINT ["vpmr"]
