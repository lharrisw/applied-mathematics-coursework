# Audio Signal Recovery — Applied Linear Algebra

An undergraduate **group project** exploring how to recover intelligible speech from an audio recording obscured by strong, rhythmic interference. The team studied the signal in both the time and frequency domains, tried several filtering methods, and documented the limitations of the resulting audio. The interference was reduced only partially; the original message was **not fully recovered**.

## My contribution

The included [`Project_2_lat.m`](Project_2_lat.m) is my individual MATLAB exploration from the project. It investigates:

- discrete Fourier transforms (FFT) and inverse transforms for inspecting and reconstructing a filtered signal;
- selection of frequency-domain components and a standard low-pass filter;
- an experimental Gaussian-like convolution filter implemented with FFTs;
- exploratory QR and singular-value decompositions applied to filtered data.

The QR/SVD portion was exploratory and is **not** presented as a validated denoising algorithm. The project supplied practical experience with numerical experimentation and evaluating approaches that did not fully solve the problem.

## Scope and reproducibility

This is **historical coursework, not production-ready signal-processing software**. The script is preserved as an individual work sample, including experimental steps and commented-out alternatives. It has not been refactored into a standalone application or independently validated for arbitrary audio. It requires an original input file named `Message.mat` containing `data` and `fs`; that file is not included, so this script is **not runnable as-is** from this repository. Some steps also require a compatible MATLAB installation and toolboxes.

The complete class report and audio recordings are intentionally **not published** here. The report was jointly authored, and including this individual script does not imply sole authorship of the overall project. Other contributors are not named on this public portfolio page.

## Context

Course: Applied Linear Algebra (undergraduate). Language: MATLAB. Subject: exploratory audio signal processing, Fourier analysis, filtering and numerical linear algebra.
