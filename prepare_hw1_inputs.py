"""Export supplied calibrated HW1 recordings as explicit HW2 inputs.

Uses only Python's standard library. Run from any directory:
    python prepare_hw1_inputs.py
No new recording, filtering, resampling, demeaning or label inference.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'hw1_inputs'
OUT.mkdir(exist_ok=True)
SOURCES = [
 ('ECG_HW1', 'biosignal_analysis/ecg_emg/results/ECG_full_channels.csv', 'ch1_mV', 'mV', .001, 1000, 144435),
 ('EMG_HW1', 'biosignal_analysis/ecg_emg/results/EMG_full_channels.csv', 'ch1_mV', 'mV', .001, 500, 18168),
 ('EEG_L03_HW1', 'biosignal_analysis/eeg/results/L03_channels_calibrated.csv', 'ch1_uV', 'uV', .000001, 200, 12910),
 ('EEG_L04_HW1', 'biosignal_analysis/eeg/results/L04_channels_calibrated.csv', 'ch1_uV', 'uV', .000001, 200, 13054),
]
manifest = {
 'course_scope': 'User-confirmed: supplied HW1 BIOPAC recordings are HW2 inputs; no additional participant acquisition is required.',
 'provenance': 'ECG/EMG electrodes connected directly to BIOPAC, as confirmed by user. HW1 identity of supplied ECG/EMG/EEG confirmed by user.',
 'operation': 'Export existing calibrated CH1 to volts. Preserve every sample and time. No filtering/resampling/demeaning. Do not apply BIOPAC calibration a second time.',
 'not_outputs': 'Prepared CSVs are HW2 input files, not measured prototype outputs, simulated outputs or proof of replay to the tablet.',
 'recording_time_limit': 'ECG/EMG contain appended runs; the stored time index is not proof of uninterrupted wall-clock recording.',
 'datasets': [],
}
for name,rel,col,unit,scale,fs,count in SOURCES:
    src=ROOT/rel
    with src.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==count
    t=[float(r['time_s']) for r in rows]
    original=[float(r[col]) for r in rows]
    volts=[v*scale for v in original]
    assert all(math.isfinite(v) for v in t+volts)
    assert all(abs(t[i]-i/fs)<1e-9 for i in range(count))
    dst=OUT/(name+'_input.csv')
    with dst.open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','input_V'])
        w.writerows((format(a,'.17g'),format(b,'.17g')) for a,b in zip(t,volts))
    with dst.open(newline='') as f:
        back=list(csv.DictReader(f))
    assert len(back)==count
    err=max(abs(float(r['input_V'])/scale-v) for r,v in zip(back,original))
    assert err<1e-10,err
    manifest['datasets'].append({
        'name':name,'source_csv':rel,'source_csv_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
        'original_calibrated_column':col,'original_unit':unit,'volts_per_original_unit':scale,
        'fs_Hz':fs,'samples':count,'length_N_over_fs_s':count/fs,'last_sample_s':t[-1],
        'output':dst.name,'output_unit':'V','output_sha256':hashlib.sha256(dst.read_bytes()).hexdigest(),
        'max_unit_roundtrip_error':err,'all_samples_preserved':True,
    })
for rel,name in [
 ('biosignal_analysis/ecg_emg/results/event_markers.csv','ECG_EMG_original_events.csv'),
 ('biosignal_analysis/eeg/results/event_markers.csv','EEG_original_events.csv'),
 ('biosignal_analysis/ecg_emg/results/channel_metadata.csv','ECG_EMG_channel_metadata.csv'),
 ('biosignal_analysis/eeg/results/channel_metadata.csv','EEG_channel_metadata.csv')]:
    shutil.copyfile(ROOT/rel,OUT/name)
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
(OUT/'README.md').write_text('''# HW1 → HW2 輸入資料

課程範圍已由使用者確認：以已提供的 HW1 BIOPAC ECG、EMG、EEG 作 HW2 輸入，不需新增人體量測。

| CSV | fs | 樣本數 | N/fs 長度 |
|---|---:|---:|---:|
| ECG_HW1_input.csv | 1000 Hz | 144435 | 144.435 s |
| EMG_HW1_input.csv | 500 Hz | 18168 | 36.336 s |
| EEG_L03_HW1_input.csv | 200 Hz | 12910 | 64.550 s |
| EEG_L04_HW1_input.csv | 200 Hz | 13054 | 65.270 s |

每個檔案只有 time_s 與 input_V 兩欄；mV 乘 0.001、uV 乘 0.000001 換為 V。既有 BIOPAC 係數已在原校準輸出使用一次，這裡不再次校準。沒有濾波、去平均、重取樣、裁段或刪點。

事件與完整通道中繼資料另附。L04 的 alpha-RMS 原單位維持在中繼資料，沒有改成 V 或頻帶功率；輸入 CSV 僅輸出 CH1。

EMG 檔實際只有臂別標記，不能自動補上 25/50/75/100%。EEG L04 A/B 未知、可能作廢檔身分未知，均維持原分析限制；不因此安排新量測。

這些是準備好的輸入，不是自製電路已量得的輸出。既有電路測試與藍牙展示保留為其原有證據。

若沿用課程的 WaveForms 檔案匯入流程，須選 input_V 欄並核對取樣率與實際電壓縮放。匯入預覽不是電壓已正確輸出的證據；不能把伏特欄當 mV，也不能把整個記錄的重播頻率當心率。此資料包不自行新增硬體量測要求。

重現：在完整原始碼根目錄執行 python prepare_hw1_inputs.py。input_manifest.json 含來源雜湊、單位換算、樣本數、時軸及往返檢查。既有分析區段與結果保持不變。
''')
print(json.dumps([{k:d[k] for k in ['name','fs_Hz','samples','length_N_over_fs_s','max_unit_roundtrip_error']} for d in manifest['datasets']],indent=2))
