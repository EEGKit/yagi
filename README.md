# yagi

Batteries-included DSP library for Rust.

## liquid-dsp, Rewritten in Rust

At the heart of this library is a full Rust implementation of [liquid-dsp](https://liquidsdr.org/), a DSP library originally written in C by Joseph Gaeddert. This library provides filters, signal modulation, mixing, error correction, and more.

liquid-dsp has more than 1000 tests, all of which are now being implemented in Yagi. The current status of the rewrite and testing efforts can be found at [LIQUID_COMPAT.md](LIQUID_COMPAT.md).

**Currently, approximately 65% of liquid-dsp's tests are passing in yagi.**

## Roadmap

Yagi is still a work in progress. The first priority is to complete the Rust rewrite of liquid-dsp. Once this is completed, we can begin investigating integrations with other libraries in order to provide SDR device access, plotting, and audio input/output. It may be useful to create some full end-to-end examples such as a stereo FM decoder.

## Created with Cursor

Development of this library shamelessly accelerated with the help of AI.
