#!/usr/bin/env python3
"""Group 2 BIOPAC EEG: fixed-window descriptive analysis, no data deletion.

Usage: python analyze_eeg.py --input raw --output results
Input can also be the supplied ZIP or its extracted directory.
BIOPAC calibrated channel.data is used ONCE; raw_data is only for rail checks.
"""
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import html
import importlib.metadata
import io
import json
import platform
import tempfile
import zipfile
from pathlib import Path

import bioread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy.signal import welch, periodogram

FS = 200.0
WELCH = dict(fs=FS, window='hann', nperseg=800, noverlap=400,
             nfft=800, detrend='constant', scaling='density', average='mean')
BANDS = {'delta': (1, 4), 'theta': (4, 8), 'alpha': (8, 13), 'beta': (13, 30)}
WINDOWS = {
    'L03': [('Closed 1', 400, 3400), ('Open', 4684, 7684), ('Closed 2', 8777, 11777)],
    'L04': [('Baseline front', 400, 3400), ('Later A (unknown)', 5054, 8054),
            ('Later B (unknown)', 9517, 12517)],
}
COLORS = ['#2367a7', '#e1812c', '#33946d']
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 180, 'figure.dpi': 110})


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with open(path, 'w', newline='', encoding='utf-8-sig') as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def integral(f, p, low, high):
    """Piecewise-linear PSD integral, inclusive exact endpoints, trapezoids.

    Boundaries here coincide with 0.25-Hz bins. Shared boundary points do not
    duplicate intervals: integral(1,4)+...+integral(13,30)==integral(1,30).
    """
    inner = (f > low) & (f < high)
    ff = np.r_[low, f[inner], high]
    pp = np.r_[np.interp(low, f, p), p[inner], np.interp(high, f, p)]
    return float(np.trapezoid(pp, ff))


def true_runs(mask):
    edges = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def longest_constant(x):
    changes = np.r_[0, np.flatnonzero(np.diff(x) != 0) + 1, len(x)]
    lengths = np.diff(changes)
    i = int(np.argmax(lengths))
    return int(changes[i]), int(lengths[i])


def quality(label, region, channel, start=0, stop=None):
    stop = len(channel.data) if stop is None else stop
    x = np.asarray(channel.data[start:stop], dtype=float)
    r = channel.raw_data[start:stop]
    is_int = np.issubdtype(r.dtype, np.integer)
    info = np.iinfo(r.dtype) if is_int else None
    rails = int(np.sum((r == info.min) | (r == info.max))) if is_int else 0
    near_rails = int(np.sum((r <= .99 * info.min) | (r >= .99 * info.max))) if is_int else 0
    pos, run = longest_constant(r)
    one_sec = [x[i:i+200] for i in range(0, len(x), 200) if len(x[i:i+200]) == 200]
    blockmeans = np.array([b.mean() for b in one_sec])
    return dict(file=label, region=region, start_s=start/FS, end_exclusive_s=stop/FS,
                samples=len(x), duration_s=len(x)/FS, nonfinite_samples=int(np.sum(~np.isfinite(x))),
                raw_min=int(r.min()), raw_max=int(r.max()),
                adc_rail_samples=rails, within_1pct_of_adc_rail_samples=near_rails,
                longest_identical_run_samples=run,
                longest_identical_run_s=run/FS, longest_identical_run_start_s=(start+pos)/FS,
                flat_1s_blocks_ptp_lt_0p1uV=sum(np.ptp(b) < .1 for b in one_sec),
                eeg_min_uV=float(x.min()), eeg_max_uV=float(x.max()),
                eeg_mean_uV=float(x.mean()), eeg_demeaned_rms_uV=float(x.std(ddof=0)),
                eeg_ptp_uV=float(np.ptp(x)), max_adjacent_step_uV=float(np.max(np.abs(np.diff(x)))),
                min_1s_mean_uV=float(blockmeans.min()), max_1s_mean_uV=float(blockmeans.max()),
                range_1s_means_uV=float(np.ptp(blockmeans)))


def channel_by_number(data, num):
    return next(c for c in data.channels if c.order_num == num)


def savefig(fig, out, name):
    # Write complete bytes atomically so consumers never see a partial image.
    for suffix in ['png', 'svg']:
        buf = io.BytesIO()
        fig.savefig(buf, format=suffix, bbox_inches='tight')
        dest = out/(name+'.'+suffix)
        tmp = out/(name+'.'+suffix+'.tmp')
        tmp.write_bytes(buf.getvalue())
        tmp.replace(dest)
    plt.close(fig)


def make_report(out, metadata, events, qrows, all_results, ratios, native_stats, versions):
    """Rebuild a Chinese report with dataset-specific, manually reviewed QC notes.

    Visual annotations below refer to the supplied Group 2 files, not to an
    automatic artifact classifier. Numerical results are regenerated from data.
    """
    blocks, md = [], []

    def heading(text, level=2):
        blocks.append(f'<h{level}>{html.escape(text)}</h{level}>')
        md.append('#'*level+' '+text+'\n')

    def para(text):
        blocks.append('<p>'+html.escape(text)+'</p>')
        md.append(text+'\n')

    def table(headers, rows):
        rows = [[str(x) for x in r] for r in rows]
        blocks.append('<div class="table-wrap"><table><thead><tr>'+''.join('<th>'+html.escape(h)+'</th>' for h in headers)+
                      '</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in r)+'</tr>' for r in rows)+
                      '</tbody></table></div>')
        md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+
                  '\n'.join('| '+' | '.join(v.replace('|','/') for v in r)+' |' for r in rows)+'\n')

    def figure(name, caption):
        data = base64.b64encode((out/(name+'.png')).read_bytes()).decode('ascii')
        blocks.append(f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>')
        md.append(f'![{caption}](results/{name}.png)\n\n{caption}\n')

    heading('Group 2 BIOPAC EEG：固定區段分析與品質審查', 1)
    para('本報告區分實際計算、目視觀察與推論。原始資料未修改；未刪除樣本、未補點、未重新選段、未額外濾波。使用者記得可能有一份作廢，但目前無法確認是哪份；本報告不替任何檔案宣告有效或作廢。')
    heading('適用範圍與品質決策')
    para('L03 有實際睜閉眼標記，但指定 Open 與 Closed 2 段有明顯可疑慢變暫態。Closed 1 亦有慢變偏移，品質未經外部驗證。因此以下三條件數字列為「固定區段、品質受限的描述性結果」，不列為已通過無偽跡審查的主要生理效應結論。此次沒有建立一組三條件皆可確認品質可接受的主要比較，也沒有為了得到預期結果另挑區段。')
    para('L04 僅保留未知事件 A、B 的時間與原始代碼，作補充品質及數值描述；不納入睜閉眼、心算或過度換氣的條件比較。後段 B 雖沒有大幅暫態，振幅小的原因與訊號來源品質仍無法僅由此檔案確認。')
    heading('檔案與原始通道')
    table(['記錄','檔名','取樣率','各通道樣本數','記錄長度'],
          [[label, next(m['filename'] for m in metadata if m['file']==label), '200 Hz',
            next(m['samples'] for m in metadata if m['file']==label),
            f"{next(m['duration_s'] for m in metadata if m['file']==label):.3f} s"] for label in ['L03','L04']])
    table(['記錄','通道','檔內名稱','原始單位欄位'],
          [[m['file'],m['channel'],m['name'],m['original_unit']] for m in metadata])
    para('BIOPAC 原生二進位檔雖無副檔名，可由 bioread 的 ACQ 解析器載入。EEG 使用已換算的 channel.data（µV），不再乘一次校正係數；raw_data 僅供整數碼與 ADC 邊界檢查。原始檔 SHA-256、比例係數、資料型別與套件版本均保存在 analysis_manifest.json / channel_metadata.csv。')
    table(['記錄','事件時間 (s)','樣本索引（0 起）','實際標記文字','代碼'],
          [[e['file'],f"{e['time_s']:.3f}",e['sample_index_zero_based'],e['original_text'] or '（空白）',e['type_code']] for e in events])
    heading('品質檢查：實際計算')
    table(['記錄／區段','時長 (s)','範圍 (µV)','ADC 邊界碼筆數','近邊界 1% 筆數','最長同值段 (s)','1秒均值最大差 (µV)'],
          [[r['file']+' / '+r['region'],f"{r['duration_s']:.3f}",f"{r['eeg_min_uV']:.3f} ～ {r['eeg_max_uV']:.3f}",
            r['adc_rail_samples'],r['within_1pct_of_adc_rail_samples'],f"{r['longest_identical_run_s']:.3f}",f"{r['range_1s_means_uV']:.3f}"] for r in qrows])
    para('兩份全程均無 NaN／Inf，未命中 int16 的精確 ADC 邊界碼；最長連續完全同值分別只有 0.010、0.015 秒。以不重疊完整 1 秒區塊之峰對峰值 < 0.1 µV 作數值平線篩查，兩份皆為 0 個；這是透明列出的輔助規則，不是經驗證的生醫品質門檻，也不能證明沒有接觸或訊號來源問題。')
    para('L04 在 64.705 秒有 1 筆原始碼 32528，接近正向 int16 邊界的 1% 範圍，換算為 198.535 µV；未形成持續平頂。此點位於所有指定區段之外。沒有命中 ADC 邊界不等於能排除上游類比放大器削波。')
    heading('品質檢查：目視觀察與位置')
    visual_rows = [
        ['L03','Closed 1：2–17 s','約 7.4–8.3 s 有慢變隆起；7.990 s 為 71.997 µV。未見長平頂或持續平線。','慢變偏移需保留限制；不是已證明乾淨的參考段'],
        ['L03','Open：23.420–38.420 s','約 23.4–24.6、28.0–32.5 s 有多個慢變暫態；23.565 s 為 −81.702 µV，32.095 s 為 132.776 µV。','有疑似偽跡；不作無偽跡主要比較'],
        ['L03','Closed 2：43.885–58.885 s','約 50.0–50.6 s 有明顯向下暫態，50.120 s 為 −135.175 µV；起點附近亦有前一暫態的回復。','有疑似偽跡；不作無偽跡主要比較'],
        ['L03','指定區段之外','約 0–1.5、43.3–43.9、63.5–64.5 s 另有大幅暫態；43.595 s 為 149.512 µV，64.190 s 為 −153.125 µV。','全程保留於圖中；不加入指定區段計算'],
        ['L04','基準前段：2–17 s','約 2–3.3、15.7–16.6 s 有尖銳或慢變暫態；2.995 s 為 −42.603 µV。','固定基準含疑似偽跡，僅補充'],
        ['L04','後段 A：25.270–40.270 s','起點附近有前一暫態的回復，之後仍有小尖波；更大的 −45.593 µV 發生在 25.205 s，早於此段。','條件未知；僅補充'],
        ['L04','後段 B：47.585–62.585 s','可見連續小幅波動，峰對峰 5.298 µV，並非持續平線；未見與末端相當的大幅暫態。','數值可描述；未知條件與小振幅來源限制仍在'],
        ['L04','指定區段之外','約 0.9、25.2、43.9、46.5、64.5–65.27 s 有暫態；末端幅度最大。','不納入指定區段；全程圖完整呈現'],
    ]
    table(['記錄','區域','目視觀察（近似範圍；極值由數值定位）','適用範圍'], visual_rows)
    write_csv(out/'quality_visual_annotations.csv', [dict(file=r[0],region=r[1],visual_observation=r[2],scope=r[3]) for r in visual_rows])
    para('推論：上述暫態的形態提示眼動、動作或電極接觸等偽跡的可能性；本資料沒有同步 EOG、動作或接觸品質紀錄，無法確診來源。這些位置是人工目視標示，不是自動刪除區間，也未用來更改固定窗。')
    figure('L03_full_EEG_events','圖 1｜L03 全程 EEG、實際事件與固定分析窗。所有樣本保留。')
    figure('L03_fixed_windows_QC','圖 2｜L03 固定三段品質圖，共用縱軸，含可疑暫態。')
    heading('固定區段與計算方法')
    table(['記錄','條件／區段','Python 半開索引','時間半開區間 (s)','樣本數'],
          [[label,n,f'[{a}:{b}]',f'[{a/FS:.3f}, {b/FS:.3f})',b-a] for label in WINDOWS for n,a,b in WINDOWS[label]])
    para('每段 15 秒、3000 samples。EEG 平均值 = mean(x)；去平均 RMS = sqrt(mean((x−mean(x))²))，等同 ddof=0 的標準差；峰對峰 = max(x)−min(x)。CH40 α RMS = sqrt(mean(CH40²))，不先去平均，另輸出去平均版本供查核。')
    para("Welch 完全採指定參數：fs=200，window='hann'，nperseg=800，noverlap=400，nfft=800，detrend='constant'，scaling='density'，average='mean'。窗長 4 秒、步距 2 秒，輸出頻率間距 0.25 Hz。每段得到 6 個重疊窗，僅作 PSD 平均，不作獨立受試者。")
    para('重要細節：3000 samples 配合上述參數，六窗涵蓋段內前 2800 samples（14 秒）；末 200 samples（1 秒）不被 Welch 使用，也不補零。這是指定參數下的演算法行為。時間域統計與 CH40 RMS 仍使用完整 3000 samples。Closed 1／Open／Closed 2 的 PSD 實際涵蓋終點分別為 16.000／37.420／57.885 秒（半開）。')
    para('頻帶功率對 PSD 作梯形積分：δ 1–4、θ 4–8、α 8–13、β 13–30 Hz，單位 µV²；PSD 單位 µV²/Hz。各端點都位於頻率格點，使用包含端點的分段線性梯形積分；相鄰頻帶只共享端點，不重複積分區間。四頻帶之和等於 1–30 Hz 積分。α 相對功率 = α 功率／1–30 Hz 功率。')
    para('未先對完整 EEG 作去趨勢或高通處理；Welch 的 constant detrend 只移除各 4 秒窗的均值，不能當作消除慢變偽跡。去平均 RMS 包含 1–30 Hz 之外的成分，因此不應期待其平方等於表列 1–30 Hz 功率。')
    heading('L03：實際計算結果（固定區段，品質受限）')
    r3 = all_results['L03']
    table(['條件','平均值 (µV)','去平均 RMS (µV)','峰對峰 (µV)','CH40 α RMS (µV)'],
          [[r['condition'],f"{r['eeg_mean_uV']:.6f}",f"{r['eeg_demeaned_rms_uV']:.6f}",f"{r['eeg_ptp_uV']:.6f}",f"{r['ch40_alpha_rms_uV']:.6f}"] for r in r3])
    table(['條件','δ (µV²)','θ (µV²)','α (µV²)','β (µV²)','1–30 Hz (µV²)','α 相對功率'],
          [[r['condition'],*[f"{r[k+'_power_uV2']:.6f}" for k in BANDS],f"{r['power_1_30Hz_uV2']:.6f}",f"{r['alpha_relative_percent']:.6f}%"] for r in r3])
    table(['比較','α 絕對功率比'],[[r['comparison'],f"{r['alpha_power_ratio']:.6f}"] for r in ratios])
    para('兩次閉眼 α 絕對功率的算術平均為 '+f"{(r3[0]['alpha_power_uV2']+r3[2]['alpha_power_uV2'])/2:.6f}"+' µV²；再除以 Open 的 α 功率。這不是把兩段閉眼串接後重新估計 PSD，也不是 α 相對功率的比值。')
    figure('L03_PSD_overlay','圖 3｜L03 三條件 PSD 疊圖；同一套指定參數，含可疑暫態，縱軸為對數。')
    figure('L03_alpha_comparison','圖 4｜L03 α 絕對及相對功率；無群體統計或推論性誤差棒。')
    para('目視觀察：PSD 在 1–4 Hz 的能量明顯，8–13 Hz 並沒有單一突出而乾淨的窄峰。實際計算顯示兩次閉眼 α 絕對功率皆高於 Open，CH40 RMS 的排序也相同，但 CH40 原廠濾波器響應與上述 8–13 Hz 積分未證實相同，僅能作輔助對照。')
    para('推論限制：數值方向與「睜眼時 α 較低」相容，但品質受限且只有單次記錄，不能僅憑功率比證實 α 阻斷或可重現的生理效應。Open 與 Closed 2 的 δ 功率高於 Closed 1，可影響 α 相對功率的分母；相對功率變化不能全數歸因於 α 本身。')
    heading('L04：補充分析，未知事件維持未知')
    figure('L04_full_EEG_events','圖 5｜L04 全程 EEG、實際標記與固定分析窗；末端大擾動在所有固定窗之外。')
    figure('L04_three_channels_events','圖 6｜L04 EEG、CH40 α、CH41 原始 α-RMS。CH41 原始單位欄位為「(uV)*秒」，定義尚待確認；圖上未作物理單位換算。')
    figure('L04_fixed_windows_QC','圖 7｜L04 固定三段 EEG；各段使用不同縱軸，以顯示小幅波動。')
    r4 = all_results['L04']
    table(['區段','平均值 (µV)','去平均 RMS (µV)','峰對峰 (µV)','CH40 α RMS (µV)'],
          [[r['condition'],f"{r['eeg_mean_uV']:.6f}",f"{r['eeg_demeaned_rms_uV']:.6f}",f"{r['eeg_ptp_uV']:.6f}",f"{r['ch40_alpha_rms_uV']:.6f}"] for r in r4])
    table(['區段','δ (µV²)','θ (µV²)','α (µV²)','β (µV²)','1–30 Hz (µV²)','α 相對功率'],
          [[r['condition'],*[f"{r[k+'_power_uV2']:.6f}" for k in BANDS],f"{r['power_1_30Hz_uV2']:.6f}",f"{r['alpha_relative_percent']:.6f}%"] for r in r4])
    para('此表的頻帶功率全部由 EEG 通道用同一 Welch 方法另算，並非讀取或換算原始 α-RMS 通道。僅描述時間區段，沒有為 A、B 賦予心算、換氣或睜眼等條件名稱。')
    table(['區段','CH41 原數值均值','CH41 原數值最小值','CH41 原數值最大值','原始單位欄位'],
          [[r['condition'],f"{r['native_channel_mean']:.6f}",f"{r['native_channel_min']:.6f}",f"{r['native_channel_max']:.6f}",r['original_unit']] for r in native_stats])
    para('CH41 只照錄為原始數值序列，定義、積分時間與原廠計算方式尚待確認；不直接改成 µV、不稱為 µV² 頻帶功率，也不與自行計算的 EEG α 功率混用。CSV 另列 native_channel_rms，僅指對此已儲存數值序列再取數學 RMS，不代表重新量得 EEG 的 RMS。')
    para('目視觀察與推論：L04 後段振幅明顯小於 L03，也小於其自身部分早段；但仍有持續波動。無同步接線、電極、受試者或增益設定紀錄，不能把低振幅逕自認定為作廢、斷線、正常閉眼或任何特定實驗效果。')
    heading('可用長度與排除規則')
    para('兩份記錄都有足夠長度容納指定的三個完整 15 秒窗，各為 45 秒、9000 samples。45 秒代表已計算資料量，不代表 45 秒全數已確認無偽跡。本次不估報一個「已驗證乾淨秒數」：Open、Closed 2、L04 基準段存在品質疑慮，其餘段也缺乏外部品質驗證。未另找看起來較符合預期的窗，未提供重新分段的主要效果估計。')
    para('最終適用範圍：L03 作有標記、但含疑似偽跡的固定區段描述；L04 作未知條件之補充與品質檢查。作廢檔案身分未確認，任何正式報告都應保留此限制，不能把此次計算包裝成已驗證有效紀錄的生理結論。')
    heading('可重現性與交付內容')
    para('執行：python -m pip install -r requirements.txt；接著 python analyze_eeg.py --input raw --output results。亦可把 --input 指向原始 ZIP 或解壓目錄。同一程式會重建數值、圖表與中文報告；人工品質註解僅在輸入雜湊與本次資料一致時產生，避免套用到其他檔案。')
    para('程式固定檢查 fs=200、每窗 3000 samples、EEG 校正關係、頻帶積分可加性，並把 Welch 與逐窗 periodogram 平均作獨立數值核對；結果存於 validation_checks.json。這些數值檢查不等同於生理訊號品質驗證。')
    para('版本：Python '+platform.python_version()+'；'+ '；'.join(k+' '+v for k,v in versions.items())+'。每段 6 個重疊窗不當作 6 名受試者；未做群體顯著性檢定、p 值、受試者層級信賴區間或推論性誤差棒。')
    table(['輸出','內容'],[
        ['L03_segment_results.csv / L03_alpha_ratios.csv','固定三條件統計、頻帶功率、CH40 RMS 與指定功率比'],
        ['L04_segment_results.csv / L04_native_alphaRMS_summary.csv','未知條件區段的 EEG 統計與原始 CH41 數值摘要'],
        ['event_markers.csv / channel_metadata.csv','原始標記文字、通道單位、樣本數、比例係數與雜湊'],
        ['PSD_all_bins.csv','兩份資料各段 0–100 Hz 的全部 PSD 頻率格點'],
        ['L03_channels_calibrated.csv / L04_channels_calibrated.csv','完整時間序列；EEG 與已知 µV 通道已校正，CH41 原始 α-RMS 保留原數值'],
        ['quality_summary.csv / quality_1s_blocks.csv','全程與固定窗的數值品質指標、不重疊1秒統計'],
        ['quality_visual_annotations.csv','人工目視疑慮位置與適用範圍'],
        ['analysis_manifest.json / validation_checks.json','參數、版本、輸入 SHA-256 與數值核對'],
        ['PNG / SVG','全程、分段品質、PSD、α 比較與 L04 三通道圖'],
        ['raw/','未修改的兩份原始檔；保留供重現，未宣告皆為有效紀錄'],
    ])
    heading('補充全程放大圖')
    figure('L03_full_record_QC_panels','補圖 1｜L03 全程分面放大，每面獨立縱軸，所有樣本皆保留。')
    figure('L04_full_record_QC_panels','補圖 2｜L04 全程分面放大，每面獨立縱軸；可看出後段不是持續平線。')
    figure('L04_PSD_overlay','補圖 3｜L04 固定區段的補充 PSD，僅為時間區段描述。')

    (out.parent/'Group2_EEG_結果與限制.md').write_text('\n'.join(md),encoding='utf-8')
    css = 'body{font-family:system-ui,"Microsoft JhengHei","Noto Sans CJK TC",sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;color:#203247;line-height:1.8}h1{font-size:29px;line-height:1.4}h2{margin-top:38px;font-size:22px;color:#165477;border-bottom:1px solid #c8dbe4;padding-bottom:8px}p{margin:14px 0}table{width:100%;border-collapse:collapse;font-size:14px;line-height:1.6}th,td{padding:10px 12px;border:1px solid #d7e1e8;text-align:left;vertical-align:top}th{background:#eaf1f6}tbody tr:nth-child(even){background:#f8fafc}.table-wrap{overflow:auto;margin:18px 0}figure{margin:26px 0}img{width:100%;height:auto;border:1px solid #dce4ea}figcaption{font-size:14px;color:#526574} @media print{body{max-width:none;margin:0;font-size:10pt}h2{break-after:avoid}figure, tr{break-inside:avoid}.table-wrap{overflow:visible}}'
    (out.parent/'Group2_EEG_Report.html').write_text('<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Group 2 BIOPAC EEG 分析</title><style>'+css+'</style><body>'+''.join(blocks)+'</body></html>', encoding='utf-8')


def event_plot(ax, label, events, xmax, short=False):
    for j, e in enumerate(events):
        # L03 has two events at zero; show the condition label and keep all in CSV.
        if label == 'L03' and e.sample_index == 0 and e.type_code == 'apnd':
            continue
        if label == 'L03':
            name = {0: 'Closed marker', 4284: 'Open marker', 8377: 'Closed marker'}.get(e.sample_index, 'Event')
        else:
            name = {0: 'Relaxed / closed marker', 4654: 'Unknown A', 9117: 'Unknown B'}.get(e.sample_index, 'Event')
        ax.axvline(e.time_index, color='#333333', ls='--', lw=.9, alpha=.8)
        if not short:
            ax.text(e.time_index + .35, .97 - .13*(j % 2),
                    name + f'\n{e.time_index:.3f} s', transform=ax.get_xaxis_transform(),
                    va='top', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', alpha=.8))
    ax.set_xlim(0, xmax)


def analyze(root, out):
    out.mkdir(parents=True, exist_ok=True)
    metadata, events, qrows, blockrows, psdrows, native_stats, checks = [], [], [], [], [], [], []
    all_results, datasets = {}, {}
    for label in ['L03', 'L04']:
        paths = [p for p in root.rglob('*') if p.is_file() and p.name.endswith('-' + label)]
        if len(paths) != 1:
            raise ValueError(f'Expected one {label} input, found {paths}')
        p = paths[0]
        d = bioread.read_file(str(p))
        datasets[label] = d
        c = channel_by_number(d, 1)
        if not all(ch.samples_per_second == FS for ch in d.channels):
            raise ValueError('Unexpected channel sample rate')
        if c.units not in ['uV', 'µV', 'μV']:
            raise ValueError(f'EEG units are not microvolts: {c.units}')
        x = np.asarray(c.data, dtype=float)  # calibrated by bioread, DO NOT scale again
        raw = c.raw_data
        if not np.allclose(x, raw.astype(float)*c.raw_scale_factor+c.raw_offset):
            raise ValueError('EEG calibration consistency check failed')
        t = np.arange(len(x)) / FS
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        for ch in d.channels:
            metadata.append(dict(file=label, filename=p.name, sha256=digest,
                                 channel=ch.order_num, name=ch.name, original_unit=ch.units,
                                 fs_Hz=ch.samples_per_second, samples=ch.point_count,
                                 duration_s=ch.point_count/ch.samples_per_second,
                                 last_sample_time_s=(ch.point_count-1)/ch.samples_per_second,
                                 storage_dtype=str(ch.dtype), raw_scale=ch.raw_scale_factor,
                                 raw_offset=ch.raw_offset))
        columns = ['time_s']+[f'ch{ch.order_num}_{"original_value" if label == "L04" and ch.order_num == 41 else "uV"}' for ch in d.channels]
        write_csv(out/f'{label}_channels_calibrated.csv',
                  (dict(zip(columns, [float(t[i])]+[float(ch.data[i]) for ch in d.channels])) for i in range(len(x))))
        for e in d.event_markers:
            events.append(dict(file=label, sample_index_zero_based=e.sample_index,
                               time_s=e.time_index, original_text=e.text, type_code=e.type_code,
                               interpretation=('unknown event A' if label == 'L04' and e.sample_index == 4654 else
                                               'unknown event B' if label == 'L04' and e.sample_index == 9117 else 'as recorded')))
        qrows.append(quality(label, 'Full record', c))
        for i in range(0, len(x), 200):
            b = x[i:i+200]
            blockrows.append(dict(file=label, start_s=i/FS, end_exclusive_s=(i+len(b))/FS,
                                  samples=len(b), mean_uV=float(b.mean()), demeaned_rms_uV=float(b.std()),
                                  min_uV=float(b.min()), max_uV=float(b.max()), ptp_uV=float(np.ptp(b))))
        rows, psds = [], []
        for name, a, b in WINDOWS[label]:
            assert b-a == 3000 and b <= len(x)
            v = x[a:b]
            f, power = welch(v, **WELCH)
            assert len(f) == 401 and np.allclose(np.diff(f), .25)
            band = {k: integral(f, power, *bounds) for k, bounds in BANDS.items()}
            total = integral(f, power, 1, 30)
            assert np.isclose(sum(band.values()), total, rtol=1e-12)
            # Independent numerical cross-check: mean of six periodograms.
            pp = [periodogram(v[j:j+800], fs=FS, window='hann', nfft=800,
                              detrend='constant', scaling='density')[1]
                  for j in range(0, len(v)-800+1, 400)]
            independent = np.mean(pp, axis=0)
            assert np.allclose(power, independent, rtol=1e-12, atol=1e-12)
            checks.append(dict(file=label, condition=name, calibration_check_passed=True,
                               all_finite=bool(np.isfinite(v).all()), samples_3000=(len(v)==3000),
                               independent_periodogram_check_passed=True,
                               max_psd_absolute_error=float(np.max(np.abs(power-independent))),
                               band_additivity_error_uV2=float(sum(band.values())-total)))
            ch40 = np.asarray(channel_by_number(d, 40).data[a:b], dtype=float)
            row = dict(file=label, condition=name, start_index=a, stop_index_exclusive=b,
                       start_s=a/FS, end_exclusive_s=b/FS, duration_s=15, samples=3000,
                       eeg_mean_uV=float(v.mean()), eeg_demeaned_rms_uV=float(np.sqrt(np.mean((v-v.mean())**2))),
                       eeg_ptp_uV=float(np.ptp(v)), **{k+'_power_uV2': val for k, val in band.items()},
                       power_1_30Hz_uV2=total, alpha_relative_fraction=band['alpha']/total,
                       alpha_relative_percent=100*band['alpha']/total,
                       ch40_alpha_rms_uV=float(np.sqrt(np.mean(ch40**2))),
                       ch40_alpha_demeaned_rms_uV=float(ch40.std()),
                       welch_windows=1+(len(v)-800)//400,
                       welch_used_samples=2800, welch_unused_trailing_samples=200,
                       welch_covered_end_exclusive_s=(a+2800)/FS,
                       scope='Fixed-window descriptive result; see quality review; not a clean confirmatory comparison')
            if label == 'L04':
                native = np.asarray(channel_by_number(d, 41).data[a:b], dtype=float)
                native_stats.append(dict(file=label, condition=name, start_s=a/FS, end_exclusive_s=b/FS,
                                         original_unit=channel_by_number(d, 41).units,
                                         native_channel_mean=float(native.mean()),
                                         native_channel_rms=float(np.sqrt(np.mean(native**2))),
                                         native_channel_min=float(native.min()), native_channel_max=float(native.max()),
                                         definition='Unconfirmed; not converted to uV or interpreted as band power'))
            rows.append(row)
            psds.append((name, f, power))
            qrows.append(quality(label, name, c, a, b))
            for hz, val in zip(f, power):
                psdrows.append(dict(file=label, condition=name, frequency_Hz=float(hz), psd_uV2_per_Hz=float(val)))
        all_results[label] = rows
        write_csv(out / f'{label}_segment_results.csv', rows)

        # Whole-record plot: full sample resolution, all values, no clipping.
        fig, ax = plt.subplots(figsize=(13, 4.7), layout='constrained')
        ax.plot(t, x, color='#273d55', lw=.6)
        for (name, a, b), color in zip(WINDOWS[label], COLORS):
            ax.axvspan(a/FS, b/FS, color=color, alpha=.13)
        event_plot(ax, label, d.event_markers, len(x)/FS)
        ax.set(xlabel='Time (s)', ylabel='EEG (µV)', title=f'{label}: full EEG, recorded markers and fixed analysis windows')
        ax.grid(alpha=.18)
        ax.legend(handles=[Patch(color=col, alpha=.25, label=n) for (n, _, _), col in zip(WINDOWS[label], COLORS)],
                  loc='lower left', ncol=3, fontsize=8)
        savefig(fig, out, f'{label}_full_EEG_events')

        # Fixed windows are retained verbatim, including suspect excursions.
        fig, axs = plt.subplots(3, 1, figsize=(13, 8), layout='constrained', sharey=(label == 'L03'))
        for ax, (name, a, b), color in zip(axs, WINDOWS[label], COLORS):
            ax.plot(t[a:b], x[a:b], color=color, lw=.7)
            ax.set(title=f'{label} / {name}: [{a}:{b}], {a/FS:.3f}–{b/FS:.3f} s', ylabel='EEG (µV)', xlim=(a/FS, b/FS))
            ax.grid(alpha=.2)
        axs[-1].set_xlabel('Time (s)')
        savefig(fig, out, f'{label}_fixed_windows_QC')

        # Zoom panels for all of the recording; independent y scales show small fluctuations.
        fig, axs = plt.subplots(7, 1, figsize=(13, 15), layout='constrained')
        for i, ax in enumerate(axs):
            a, b = i*2000, min((i+1)*2000, len(x))
            ax.plot(t[a:b], x[a:b], lw=.65, color='#273d55')
            ax.set(ylabel='µV', title=f'{label}: {a/FS:.1f}–{b/FS:.3f} s', xlim=(a/FS, b/FS))
            ax.grid(alpha=.2)
        axs[-1].set_xlabel('Time (s); panels use independent y scales')
        savefig(fig, out, f'{label}_full_record_QC_panels')

        fig, ax = plt.subplots(figsize=(9, 5), layout='constrained')
        for (name, f, power), color in zip(psds, COLORS):
            ax.semilogy(f, power, color=color, lw=1.6, label=name)
        for low, high in BANDS.values():
            ax.axvline(low, color='gray', alpha=.3, lw=.7)
        ax.axvline(30, color='gray', alpha=.3, lw=.7)
        visible = np.concatenate([power[(f>=1)&(f<=30)] for _, f, power in psds])
        ax.set(xlim=(1, 30), ylim=(visible.min()*.65, visible.max()*1.5), xlabel='Frequency (Hz)', ylabel='PSD (µV²/Hz)',
               title=f'{label}: Welch EEG PSD / fixed 15-s windows / no artifact removal')
        ax.legend()
        ax.grid(alpha=.2, which='both')
        savefig(fig, out, f'{label}_PSD_overlay')

    # L03 ratios use arithmetic average of the two CLOSED powers.
    rows = all_results['L03']
    ap = np.array([r['alpha_power_uV2'] for r in rows])
    ratios = [dict(comparison='Closed 1 / Open', alpha_power_ratio=ap[0]/ap[1]),
              dict(comparison='Closed 2 / Open', alpha_power_ratio=ap[2]/ap[1]),
              dict(comparison='Mean(Closed 1, Closed 2) / Open', alpha_power_ratio=((ap[0]+ap[2])/2)/ap[1])]
    write_csv(out/'L03_alpha_ratios.csv', ratios)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6), layout='constrained')
    for ax, key, title, unit in [(axs[0], 'alpha_power_uV2', 'Absolute alpha power (8–13 Hz)', 'Power (µV²)'),
                                 (axs[1], 'alpha_relative_percent', 'Relative alpha power / 1–30 Hz', 'Relative power (%)')]:
        vals = [r[key] for r in rows]
        ax.bar([r['condition'] for r in rows], vals, color=COLORS, width=.62)
        for i, val in enumerate(vals):
            ax.text(i, val, f'{val:.3f}', ha='center', va='bottom')
        ax.set(ylabel=unit, title=title, ylim=(0, max(vals)*1.18))
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    fig.suptitle('L03: fixed windows include suspect transients; descriptive results')
    savefig(fig, out, 'L03_alpha_comparison')

    d = datasets['L04']
    fig, axs = plt.subplots(3, 1, figsize=(13, 9), sharex=True, layout='constrained')
    for ax, num, title, unit in zip(axs, [1, 40, 41], ['EEG', 'CH40 alpha', 'CH41 original alpha-RMS'],
                                  ['µV', 'µV', 'Original CH41 values\n(unit field in report / CSV)']):
        ch = channel_by_number(d, num)
        tt = np.arange(len(ch.data))/FS
        ax.plot(tt, ch.data, color='#273d55', lw=.65)
        for (name, a, b), color in zip(WINDOWS['L04'], COLORS):
            ax.axvspan(a/FS, b/FS, color=color, alpha=.13)
        event_plot(ax, 'L04', d.event_markers, len(ch.data)/FS, short=(num != 1))
        ax.set(ylabel=unit, title=title)
        ax.grid(alpha=.2)
    axs[-1].set_xlabel('Time (s)')
    fig.suptitle('L04: unknown events remain unknown; CH41 definition unconfirmed')
    savefig(fig, out, 'L04_three_channels_events')

    write_csv(out/'channel_metadata.csv', metadata)
    write_csv(out/'event_markers.csv', events)
    write_csv(out/'quality_summary.csv', qrows)
    write_csv(out/'quality_1s_blocks.csv', blockrows)
    write_csv(out/'PSD_all_bins.csv', psdrows)
    write_csv(out/'L04_native_alphaRMS_summary.csv', native_stats)
    versions = {key: importlib.metadata.version(key) for key in ['numpy','scipy','matplotlib','bioread']}
    manifest = dict(python=platform.python_version(), packages=versions, fs_Hz=FS,
                    welch=WELCH, bands_Hz=BANDS, windows=WINDOWS,
                    band_integration='piecewise-linear trapezoids including exact band endpoints',
                    eeg_scaling='bioread channel.data, already calibrated to uV; no additional scaling',
                    no_artifact_removal=True, no_hypothesis_testing=True,
                    discarded_record_identity='unknown', input_channels=metadata,
                    l03_ratios=ratios)
    (out/'analysis_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
    (out/'validation_checks.json').write_text(json.dumps(checks, indent=2, ensure_ascii=False), encoding='utf-8')
    known = {'L03': 'd6b455c42473900460ecb5b06019312c68820cff4497b3fe88034f5ef18ef4bd',
             'L04': '905a5cc835ddb186d6364cfeec079dcffe701af8c54cf6a8b298a98225d15a55'}
    if all(next(m['sha256'] for m in metadata if m['file']==label) == digest for label, digest in known.items()):
        make_report(out, metadata, events, qrows, all_results, ratios, native_stats, versions)
    else:
        print('Input hashes differ: numerical results produced; dataset-specific manual QC report not applied.')
    print(json.dumps(dict(L03=all_results['L03'], L04=all_results['L04'], ratios=ratios), ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(__file__).parent/'raw')
    parser.add_argument('--output', type=Path, default=Path(__file__).parent/'results')
    args = parser.parse_args()
    if args.input.is_file() and zipfile.is_zipfile(args.input):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with zipfile.ZipFile(args.input) as z:
                for member in z.infolist():
                    if not (root/member.filename).resolve().is_relative_to(root.resolve()):
                        raise ValueError('Unsafe ZIP path')
                z.extractall(root)
            analyze(root, args.output)
    else:
        analyze(args.input, args.output)


if __name__ == '__main__':
    main()
