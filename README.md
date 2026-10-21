# Group 2 — Biomedical Experiment 2

AD620 / LM358 analog signal chain, Arduino Uno ADC, HMSoft BLE and a custom MIT App Inventor waveform viewer; reproducible offline HW1 BIOPAC EEG / ECG / EMG analysis.

## Course code entry points

| Item | Source |
|---|---|
| Arduino ADC / BLE | [Group2_Biosignal.ino](arduino/Group2_Biosignal/Group2_Biosignal.ino) |
| EEG processing and four-band comparison | [analyze_eeg.py](biosignal_analysis/eeg/analyze_eeg.py) |
| Custom Android App project | [Group2_Biosignal.aia](app_development/Group2_Biosignal.aia) |
| App packet parser and plot | [waveform.html](app_development/waveform.html) |
| ECG / EMG analysis | [analyze_ecg_emg.py](biosignal_analysis/ecg_emg/analyze_ecg_emg.py) |
| HW1 input unit conversion | [prepare_hw1_inputs.py](prepare_hw1_inputs.py) |

## Hardware / firmware

- Analog chain: AD620 → U1 RC bandpass → passive Twin-T → U2 buffer / non-inverting postamp → U3 positive averager → Uno A0.
- AD620 RG 470 Ω (nominal gain ≈106.11). U2-B Rf ≈3.4 kΩ, Rg ≈910 Ω (gain ≈4.736). U3 uses two 2.46 kΩ branches, nominal VA0 ≈2.5 + 0.5 Vpost.
- HMSoft carrier TXD → Uno D10. Uno D11 → 470 Ω → J → RXD; J → 390 Ω + 470 Ω series → GND. The analog node S is separate from J.
- Uno USB supply; Uno GND directly to AD3 GND. Analog rails nominal ±5 V; Uno 5 V is separate from AD3 V+.
- USB Serial 115200 baud; BLE UART 9600 baud. Frames: `*Q<raw ADC integer>*`; nominal sampling interval 10 ms. Actual sampling rate / jitter / loss are not calibrated.

This is a generator-driven bench setup. Follow the documented laboratory workflow; it is not a validated direct human acquisition interface. Existing human recordings were acquired using BIOPAC in HW1 and are reused as HW2 inputs. No additional participant acquisition is required.

## Run the Android App

Import `app_development/Group2_Biosignal.aia` into MIT App Inventor; retain the bundled BLE extension (20240822). Rename an existing same-name project first. Connect via AI Companion, scan, select HMSoft, then verify the waveform. Build an APK using MIT App Inventor's Build menu when ready.

The app supports full streaming frame parsing across BLE notifications, the latest 100 ADC samples, minimum/maximum, pause/resume, clear, raw diagnostics and CSV export/sharing. X is received-sample index, Y is 0–1023 ADC counts. The cumulative valid-frame count is not a packet-loss metric. The successful Companion bench display showed a 400–645-count window; this is not a calibrated voltage or physiological output. Standalone APK and Android sharing need their own confirmation.

The `.aia` files use a generic `ai_group2` namespace for publication. The bundled extension binaries are unchanged. No APK is included.

## Reproduce analysis locally

Python 3.12 is the recorded environment. Install each analysis folder's pinned `requirements.txt` in a virtual environment. Raw BIOPAC recordings are intentionally not distributed here; place your authorized HW1 files under the corresponding `raw/` directory.

```bash
python -m pip install -r biosignal_analysis/eeg/requirements.txt
python biosignal_analysis/eeg/analyze_eeg.py --input biosignal_analysis/eeg/raw --output biosignal_analysis/eeg/results
python biosignal_analysis/ecg_emg/analyze_ecg_emg.py --input biosignal_analysis/ecg_emg/raw --output biosignal_analysis/ecg_emg/results
python prepare_hw1_inputs.py
```

Input filenames may be extensionless BIOPAC ACQ files. EEG analysis locates the L03/L04 recordings; ECG and EMG code identifies the actual channels/recordings rather than assuming that a filename containing ECG is ECG. `prepare_hw1_inputs.py` runs after both analyses and uses their calibrated full-channel CSVs. It converts mV / µV to V once, retaining all samples and timestamps.

EEG: fs=200 Hz, Hann Welch windows 800, overlap 400, FFT 800, constant detrend, density scaling, mean averaging. Bands δ 1–4, θ 4–8, α 8–13, β 13–30 Hz; α relative power uses 1–30 Hz as denominator. L03 fixed slices: [400:3400], [4684:7684], [8777:11777]. L04 later events remain unknown. No artifact deletion or group significance test is added. EMG has arm labels, not verified individual force-level labels.

## Rebuild / test App code

Run from this repository root:

```bash
python app_development/build_project.py
node app_development/test_waveform.js
python app_development/test_native_blocks.py
```

The builder takes `upload/Group2_Biosignal.aia` (the minimal, sanitized original project) and produces the complete `.aia`. Tests cover 100 random fragmentation runs of 1000 samples each and native-block control flow with mocks. These tests do not emulate Android BLE, MIT compilation or layout rendering.

## Data and evidence boundaries

- No HW1 participant recordings, sample-level physiological CSVs, course handouts, full report PDF or tablet screenshots are distributed in this code repository.
- Input data, bench stage measurements and App observations remain distinct evidence. No synthetic sample is presented as a measurement.
- The identity of a possibly discarded EEG recording is unknown; no record is discarded because its alpha trend differs from expectation.
- GitHub publication does not resolve existing acquisition / artifact / bandwidth limitations.

See [CONTRIBUTIONS.md](CONTRIBUTIONS.md) and [NOTICE.md](NOTICE.md).
