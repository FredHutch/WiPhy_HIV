FROM ubuntu:20.04

# Install pre-requisites
RUN apt-get update
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y \
    build-essential \
    libgtkglext1-dev \
    libgsl-dev \
    libglu-dev \
    libgtk-3-dev \
    libpango1.0-dev \
    libpangocairo-1.0-0 \
    libgtk2.0-dev \
    libgtk2.0-0 \
    libgtk-3-0


ADD . /usr/local/src
WORKDIR /usr/local/src/Linux
RUN export CPATH=$CPATH:/usr/include/gtk-2.0 \
    export CPATH=$CPATH:/usr/include/glib-2.0/ \
    export CPATH=$CPATH:/usr/lib/glib-2.0/include/ \
    export CPATH=$CPATH:/usr/include/pango-1.0/ \
    export CPATH=$CPATH:/usr/lib/gtk-2.0/include/ \
    make
