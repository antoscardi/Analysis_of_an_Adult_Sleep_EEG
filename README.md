# 💤 Analysis of an Adult Sleep EEG

## 📌 Overview
This project analyzes overnight EEG recordings 🧠 to identify sleep stages and patterns. The MATLAB script applies various signal processing techniques 🎛️ to extract meaningful brain activity trends over a night of sleep.

## 🔍 Features
- **Preprocessing 🎚️:** Down-sampling, mean subtraction, and anti-aliasing filtering.
- **Noise Cancellation 🧹:** Band-pass filtering to isolate EEG frequency bands.
- **Wave Extraction 🌊:** Band-pass filters for Alpha, Beta, Delta, and Theta waves.
- **Sleep Phase Analysis 🛏️:** Spectrogram visualization and estimation of sleep stages.
- **Hypnogram Generation 📊:** A speculative sleep stage representation based on spectral analysis.

## 🚀 How to Run
### Prerequisites
- MATLAB 💻
- Signal Processing Toolbox 🛠️

### Steps
1. Load the EEG dataset 📂.
2. Run the MATLAB script ▶️.
3. Follow the steps:
   - Down-sample the signal ⏬.
   - Apply band-pass filtering for noise removal 🎛️.
   - Extract characteristic brain waves 🧠.
   - Generate and analyze spectrograms 📊.
4. The script will output:
   - Processed EEG signals 📉.
   - Spectrograms for different brain waves 🔬.
   - A speculative hypnogram illustrating sleep stage transitions 💤.

## 📑 Output
- **Filtered EEG Signals 🎚️:** Displays cleaned EEG signals.
- **Spectrograms 🔬:** Shows time-frequency representation of EEG waves.
- **Hypnogram 🏷️:** A speculative representation of sleep stages based on EEG analysis.

## 📖 Methodology
1. **Preprocessing 🧼:**
   - Mean subtraction for signal centering.
   - Anti-aliasing filtering to remove unwanted frequencies.
   - Down-sampling to reduce computational complexity.
2. **Noise Cancellation 🧹:**
   - Band-pass filtering to remove baseline wander and high-frequency noise.
3. **Wave Extraction 🌊:**
   - Extracting Alpha (8-13 Hz), Beta (14-35 Hz), Delta (0.5-4 Hz), and Theta (4-7 Hz) waves.
   - Using band-pass filters to isolate each wave type.
4. **Spectrogram Analysis 📊:**
   - Visualizing brain activity over time.
   - Identifying transitions between sleep phases.
5. **Hypnogram Generation 💤:**
   - Constructing a speculative sleep stage timeline.

## 🎯 Example Usage
Upon running the script ▶️, the user follows on-screen prompts to:
- Filter and analyze EEG data 🧠.
- Identify sleep stage transitions 📉.
- Generate an overview of the night's sleep phases 💤.

## 👥 Authors
- Guglielmo Bruno (Politecnico di Milano) 🎓
- Mirko Coggi (Politecnico di Milano) 🎓
- Antonio Scardino (Politecnico di Milano) 🎓

For any inquiries, contact:
- 📧 guglielmo.bruno@mail.polimi.it
- 📧 mirko.coggi@mail.polimi.it
- 📧 antonio.scardino@mail.polimi.it
