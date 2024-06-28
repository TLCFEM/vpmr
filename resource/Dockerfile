FROM debian AS builder
RUN apt-get update && apt-get install -y cmake git g++ libtbb-dev libmpfr-dev libgmp-dev
RUN git clone --depth 1 https://github.com/TLCFEM/vpmr.git
WORKDIR /vpmr
RUN git submodule update --init --recursive
RUN sed -i 's/3.24/3.18/g' CMakeLists.txt
RUN cd eigen && git apply --ignore-space-change --ignore-whitespace ../patch_size.patch && cd ..
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

FROM debian AS runtime
RUN apt-get update && apt-get install -y libtbb12 libmpfr6 libgmp10
COPY --from=builder /vpmr/build/vpmr /usr/bin/vpmr
ENTRYPOINT ["/usr/bin/vpmr"]
CMD ["-h"]
