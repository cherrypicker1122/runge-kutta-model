FROM gcc:4.9
COPY . /usr/src/rk4
WORKDIR /usr/src/rk4
RUN /usr/bin/g++-12 -o rk4 twobody.cpp
CMD ["./rk4"]
