FROM alpine:3 AS builder

RUN apk add make cmake git g++ onetbb-dev mpfr-dev gmp-dev

RUN git clone --depth 1 https://github.com/TLCFEM/vpmr.git
WORKDIR /vpmr
RUN git submodule update --init --recursive && cd eigen && git apply --ignore-space-change --ignore-whitespace ../patch_size.patch && cd .. && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

FROM alpine:3 AS runtime
RUN apk add onetbb mpfr
COPY --from=builder /vpmr/build/vpmr /usr/bin/vpmr

ENTRYPOINT ["/usr/bin/vpmr"]
CMD ["-h"]
