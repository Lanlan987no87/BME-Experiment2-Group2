#!/usr/bin/env python3
"""Group 2 ECG / EMG: read BIOPAC channel metadata, use recorded conditions.

python analyze_ecg_emg.py --input raw --output results
Requires Python 3.12, numpy, scipy, matplotlib, bioread.
No clinical diagnosis; no unlabeled force/load conditions are invented.
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
from pathlib import Path

import bioread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfiltfilt, find_peaks

COLORS = ['#246ca8', '#299777', '#a17ab5', '#df8743']
ECG_LABELS = ['Supine', 'Seated', 'Deep breathing', 'Post exercise']
ECG_ZH = ['仰臥', '坐姿', '深呼吸', '運動後']
EMG_LABELS = ['Dominant arm', 'Nondominant arm']
EMG_ZH = ['主力臂', '非主力臂']
DETECTOR = dict(bandpass_Hz=[8,30], butter_order=3, min_peak_distance_s=.300,
                prominence_mV=.12, raw_peak_refine_halfwidth_s=.040,
                analysis_boundary_guard_s=1.0)
KNOWN = {
 'ECG': '69d472e539a7caf6d013adb4cd899e0976eb8218888b78d3d74f1706900a6931',
 'EMG': '43c3287a290e9a4de3b4bea94bb2ec2abc00d18b51ba5fb27fdc293e5452adf4',
}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':170,'figure.dpi':110})


def write_csv(path, rows):
    rows=list(rows)
    if not rows:return
    with open(path,'w',encoding='utf-8-sig',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def savefig(fig,out,name):
    for fmt in ['png','svg']:
        b=io.BytesIO();fig.savefig(b,format=fmt,bbox_inches='tight')
        final=out/(name+'.'+fmt);tmp=out/(name+'.'+fmt+'.tmp')
        tmp.write_bytes(b.getvalue());tmp.replace(final)
    plt.close(fig)


def qc(label,condition,ch,start,stop):
    x=np.asarray(ch.data[start:stop],dtype=float);r=ch.raw_data[start:stop];fs=ch.samples_per_second
    info=np.iinfo(r.dtype)
    change=np.r_[0,np.flatnonzero(np.diff(r)!=0)+1,len(r)]
    run=int(np.diff(change).max())
    n=int(fs);blocks=[x[j:j+n] for j in range(0,len(x)-n+1,n)]
    means=np.array([v.mean() for v in blocks])
    return dict(file=label,condition=condition,start_index=start,stop_index_exclusive=stop,
                start_s=start/fs,end_exclusive_s=stop/fs,samples=len(x),duration_s=len(x)/fs,
                original_unit=ch.units,nonfinite_samples=int((~np.isfinite(x)).sum()),
                adc_rail_samples=int(((r==info.min)|(r==info.max)).sum()),
                within_1pct_of_adc_rail_samples=int(((r<=.99*info.min)|(r>=.99*info.max)).sum()),
                longest_exact_constant_run_s=run/fs,
                flat_1s_blocks_ptp_lt_one_raw_lsb=int(sum(np.ptp(b)<abs(ch.raw_scale_factor) for b in blocks)),
                raw_min=int(r.min()),raw_max=int(r.max()),min_mV=float(x.min()),max_mV=float(x.max()),
                mean_mV=float(x.mean()),demeaned_rms_mV=float(x.std(ddof=0)),ptp_mV=float(np.ptp(x)),
                range_1s_means_mV=float(np.ptp(means)),max_adjacent_step_mV=float(np.abs(np.diff(x)).max()))


def r_peaks(raw,fs,prominence):
    """Dataset-specific positive R detector; no RR interval deletion.

    Filter separately within each appended recording to avoid filtering across
    concatenation boundaries. The calibrated RAW local maximum is the R time.
    """
    sos=butter(DETECTOR['butter_order'],DETECTOR['bandpass_Hz'],
               btype='bandpass',fs=fs,output='sos')
    filtered=sosfiltfilt(sos,raw)
    peaks,_=find_peaks(filtered,distance=round(DETECTOR['min_peak_distance_s']*fs),
                       prominence=prominence)
    half=round(DETECTOR['raw_peak_refine_halfwidth_s']*fs)
    refined=[]
    for peak in peaks:
        a=max(0,int(peak)-half);b=min(len(raw),int(peak)+half+1)
        refined.append(a+int(np.argmax(raw[a:b])))
    return np.unique(refined).astype(int),filtered


def decorate_events(ax,events,fs,duration,labels,annotate=True):
    for i,e in enumerate(events):
        ax.axvline(e.time_index,color='#444444',lw=.8,ls='--')
        if annotate:
            ax.text(e.time_index+.5,.97-.15*(i%2),labels[i]+f'\n{e.time_index:.3f} s',
                    transform=ax.get_xaxis_transform(),va='top',fontsize=8,
                    bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    ax.set_xlim(0,duration);ax.grid(alpha=.18)


def make_report(out,metadata,events,quality,ecg,emg,ratios,breath,sensitivity,versions):
    blocks=[];md=[]
    def h(text,level=2):
        blocks.append(f'<h{level}>{html.escape(text)}</h{level}>');md.append('#'*level+' '+text+'\n')
    def p(text):
        blocks.append('<p>'+html.escape(text)+'</p>');md.append(text+'\n')
    def table(headers,rows):
        rows=[[str(v) for v in r] for r in rows]
        blocks.append('<div class="table-wrap"><table><thead><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in headers)+
                      '</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(c)+'</td>' for c in r)+'</tr>' for r in rows)+
                      '</tbody></table></div>')
        md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+
                  '\n'.join('| '+' | '.join(c.replace('|','/') for c in r)+' |' for r in rows)+'\n')
    def fig(name,caption):
        b=base64.b64encode((out/(name+'.png')).read_bytes()).decode('ascii')
        blocks.append(f'<figure><img src="data:image/png;base64,{b}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>')
        md.append(f'![{caption}](results/{name}.png)\n\n{caption}\n')

    h('Group 2 BIOPAC ECG／EMG：實際條件分析與品質核對',1)
    p('分析範圍以檔內通道與實際標記為準。報告分開列出實際計算、目視觀察與推論。ECG 與 EMG 都是 BIOPAC／BSL 原生記錄。使用者於 2026-10-20 確認電極直接接 BIOPAC，未經自製 AD620／LM358 前端；此為使用者確認的量測來源。使用者進一步確認這些是 HW1 資料，作 HW2 輸入，不需額外人體量測。')
    h('先確認檔案身分、通道與單位')
    table(['檔名','實際類型','通道','檔內名稱','原始單位','取樣率','樣本數','時長'],
          [[m['filename'],m['file'],m['channel'],m['name'],m['original_unit'],f"{m['fs_Hz']:g} Hz",m['samples'],f"{m['duration_s']:.3f} s"] for m in metadata])
    p('重要辨識：1151_g02_ecg02-L01 雖然檔名有 ecg，實際 CH1 名稱是 EMG，內嵌文件也為 L01 肌電圖。本次把它作 EMG 分析，沒有當成第二份 ECG。CH40 積分 EMG 的原始單位 mV-sec 保留，不改稱原始 EMG RMS。')
    p('校正：主通道使用 bioread 的已校正 channel.data（mV），未再乘一次比例係數；並核對 channel.data = raw_data × scale + offset。原始檔、SHA-256、校正係數與完整通道序列均附於交付檔。')
    h('實際條件與事件標記')
    manual=[e for e in events if e['type_code']!='nrto']
    table(['資料','檔內時間 (s)','樣本索引（0起）','標記文字','代碼'],
          [[e['file'],f"{e['time_s']:.3f}",e['sample_index_zero_based'],e['original_text'],e['type_code']] for e in manual])
    p('ECG 另有 229 個 nrto 自動「新次率」標記，全部原樣輸出至 event_markers.csv；它們是軟體衍生輸出，不是 229 次人工條件切換，也不當作 R 波真值。檔內時間是多段附加記錄的串接時間，不能據此認定姿勢改變、運動或休息的過程皆連續記錄。')
    p('ECG 量測狀態可辨識為仰臥、坐姿、深呼吸、運動後，並有 5 組開始吸氣／開始呼氣標記。EMG 只有主力臂與非主力臂；沒有每次收縮、放鬆、重量、負載或左右手標記。主力臂不逕自等同右手，非主力臂不逕自等同左手。')
    h('波形品質：實際數值檢查')
    table(['資料／區段','範圍 (mV)','ADC邊界碼筆數','近ADC邊界1%筆數','最長連續同值 (s)','1秒均值範圍差 (mV)'],
          [[q['file']+' / '+q['condition'],f"{q['min_mV']:.6f} ～ {q['max_mV']:.6f}",q['adc_rail_samples'],q['within_1pct_of_adc_rail_samples'],
            f"{q['longest_exact_constant_run_s']:.3f}",f"{q['range_1s_means_mV']:.6f}"] for q in quality])
    p('兩份主通道皆無 NaN／Inf，沒有命中或進入 int16 ADC 邊界的 1% 範圍。全程最長連續完全同值為 ECG 0.012 秒、EMG 0.008 秒，未見持續平線。另以完整不重疊 1 秒的峰對峰值小於一個原始 ADC LSB 作數值平線篩查；結果見 quality_summary.csv。未命中 ADC 邊界不能排除上游類比處理的所有限制。')
    h('波形品質：目視觀察與適用範圍')
    annotations=[
        ['ECG','仰臥','R 波尖銳且可分辨，後方有較寬的 T 波；原廠心率通道在此段明顯高於由 R–R 重算的值。','採重新定位的 R 波計算，不使用原廠 CH40 作主要心率'],
        ['ECG','坐姿／深呼吸／運動後','約 29.5–30.8、42–44、49–50、63–64、90–91、94–95、99–100、117–123 秒等處有基線位移或暫態；部分深呼吸段可見週期性基線起伏。','R 峰仍可追蹤；適用於節拍時間分析，不延伸至 ST、P／T 振幅或精細形態分析'],
        ['EMG','主力臂','可見約 2–4、6.3–7.6、10–11.7、14–15.8 秒的活動群及間隔。','僅作目視活動分布描述；未用這些近似時間切成正式收縮條件'],
        ['EMG','非主力臂','可見約 19.5–22、24–26、28–30、32–34 秒的活動群；約 31 秒後低活動間隔有慢變偏移，最後一群含大幅尖波。','完整臂別區段照算；峰對峰對單一尖波敏感，不等同肌力'],
    ]
    table(['資料','區域','目視觀察（近似時間）','適用範圍'],annotations)
    write_csv(out/'quality_visual_annotations.csv',[dict(file=r[0],region=r[1],visual_observation=r[2],scope=r[3]) for r in annotations])
    p('推論：基線位移可能與呼吸、動作或接觸變化相關，但沒有同步影像、呼吸帶或接觸紀錄，不能確定來源。EMG 的活動群不能單靠波形重新命名為某個指定負載或最大握力，也不把高振幅本身當成偽跡或生理強度的充分證據。')
    fig('ECG_full_events_quality','圖 1｜完整 ECG、條件標記、原廠與重算心率，以及不重疊1秒均值。所有原始樣本保留。')
    h('ECG：R 波定位與分析區段')
    p('各附加記錄各自處理，避免在串接邊界跨段濾波。R 波偵測使用 scipy.signal.butter 的 N=3、8–30 Hz bandpass 設定，再以 sosfiltfilt 前後向零相位處理；N=3 是設計參數，bandpass 的單向總階數為6。在正向濾波訊號上用 find_peaks，最小距離 0.300 秒、prominence 0.12 mV。每個候選點再回到原始已校正 ECG，於 ±0.040 秒找局部最大值，所得原始樣本時間才作 R 波時間。這是針對本次正向 R 波資料的設定，不是通用或經臨床驗證的偵測器。')
    p('為降低附加記錄起訖與濾波邊界效應，四個實際條件皆固定排除首尾各 1 秒；規則套用於所有條件，不依心率高低改窗。R–R 只有當兩端 R 波都落在同一分析區段時才納入，不跨條件接續。所有偵測點（含邊界未納入點）均保留於 ECG_R_peaks.csv；沒有人工增刪 R 點或刪除 R–R。')
    table(['條件','原記錄半開時間 (s)','分析 Python 半開索引','分析半開時間 (s)','R 點數','R–R 數'],
          [[r['condition'],f"[{r['record_start_s']:.3f}, {r['record_end_exclusive_s']:.3f})",
            f"[{r['analysis_start_index']}:{r['analysis_stop_index_exclusive']}]",
            f"[{r['analysis_start_s']:.3f}, {r['analysis_end_exclusive_s']:.3f})",r['R_peaks'],r['RR_intervals']] for r in ecg])
    p('數值驗證：prominence 0.12 與 0.15 mV 在四個分析區段得到完全相同的 R 波時間；0.10 mV 在深呼吸段產生額外候選點，故不採用。另核對 R–R 數 = R 點數−1，R–R 加總 = 最末 R 時間−第一 R 時間。全程 R 點疊圖已目視檢查；這些檢查支持本次節拍定位，但不是外部人工真值驗證。')
    fig('ECG_R_peak_examples','圖 2｜四條件原始 ECG 局部與 R 峰定位。原始波形不因偵測濾波而被取代。')
    h('ECG：實際計算結果')
    p('R–R = 相鄰 R 波樣本索引差／fs；主要心率 = 60／平均 R–R（秒）。另輸出各間期瞬時心率 60／R–R 與其算術平均。兩種心率平均定義不同，報告主表明確採前者。R–R 標準差使用 ddof=1，只描述本段間期分散，不將它宣稱為經品質控制的 NN 或正式 HRV 指標。')
    table(['條件','平均 R–R (ms)','R–R 標準差 (ms)','R–R 最小～最大 (ms)','60／平均R–R (bpm)','瞬時HR平均 (bpm)'],
          [[r['condition'],f"{r['RR_mean_ms']:.6f}",f"{r['RR_std_ms_ddof1']:.6f}",f"{r['RR_min_ms']:.3f} ～ {r['RR_max_ms']:.3f}",
            f"{r['HR_from_mean_RR_bpm']:.6f}",f"{r['mean_instantaneous_HR_bpm']:.6f}"] for r in ecg])
    fig('ECG_RR_HR_timeline','圖 3｜R–R 與瞬時心率時間圖。各條件分別連線，避免跨段串接；深呼吸區的細直線為實際呼吸相位標記。')
    fig('ECG_condition_comparison','圖 4｜同一記錄中四個實際條件的心率／平均 R–R 描述性比較，沒有群體推論誤差棒。')
    h('ECG：為何不直接採原廠 CH40 心率')
    table(['條件','原廠CH40區段均值 (bpm)','由R–R重算主要心率 (bpm)'],
          [[r['condition'],f"{r['native_CH40_mean_bpm']:.6f}",f"{r['HR_from_mean_RR_bpm']:.6f}"] for r in ecg])
    p('實際計算：仰臥段原廠 CH40 均值為 '+f"{ecg[0]['native_CH40_mean_bpm']:.3f}"+' bpm，由原始 R 波間隔重算為 '+f"{ecg[0]['HR_from_mean_RR_bpm']:.3f}"+' bpm。目視相鄰 R 波約相隔 0.8 秒，但原廠通道多在約 120–180 bpm 間跳動。新次率標記只代表自動輸出更新時間，不等於已確認的 QRS 觸發時間。推論：原廠設定可能偵測到非 R 成分或有其他設定問題；沒有內部設定，不能僅憑新次率時間證實 T 波誤判。此次採重算 R–R，保留原廠通道作查核。')
    fig('ECG_native_HR_discrepancy','圖 5｜仰臥 2–6 秒原始 ECG 與心率對照。上圖橘色直線為原廠新次率標記，紅點是重算 R 峰。')
    h('ECG：實際呼吸標記的補充摘要')
    table(['相位序號','實際標記','標記時間 (s)','分析相位時間 (s)','完整R–R數','平均R–R (ms)','心率 (bpm)'],
          [[r['phase_number'],r['recorded_phase_marker'],f"{r['marker_time_s']:.3f}",f"[{r['analysis_phase_start_s']:.3f}, {r['analysis_phase_end_exclusive_s']:.3f})",r['RR_intervals'],
            f"{r['RR_mean_ms']:.3f}" if r['RR_mean_ms']!='' else '無完整間期',
            f"{r['HR_from_mean_RR_bpm']:.3f}" if r['HR_from_mean_RR_bpm']!='' else '不計算'] for r in breath])
    p('各相位僅取兩個 R 點皆在該相位內的完整間期，不把跨相位間期硬分配給吸氣或呼氣。少量間期及標記時差會影響相位摘要；不據此宣稱呼吸相位效果具有統計顯著性。')
    h('EMG：僅比較可辨識的臂別條件')
    p('計算使用原始、已校正的 EMG CH1（mV），不額外濾波或整流。去平均 RMS = sqrt(mean((x−mean(x))²))；峰對峰 = max(x)−min(x)。條件邊界只採實際主力臂／非主力臂標記，完整區段所有樣本皆保留。')
    table(['實際條件','Python半開索引','半開時間 (s)','樣本數','時長 (s)','平均值 (mV)','去平均RMS (mV)','峰對峰 (mV)'],
          [[r['condition'],f"[{r['start_index']}:{r['stop_index_exclusive']}]",f"[{r['start_s']:.3f}, {r['end_exclusive_s']:.3f})",r['samples'],f"{r['duration_s']:.3f}",
            f"{r['mean_mV']:.6f}",f"{r['demeaned_rms_mV']:.6f}",f"{r['ptp_mV']:.6f}"] for r in emg])
    p('非主力臂／主力臂：去平均 RMS 比為 '+f"{ratios[0]['demeaned_rms_ratio']:.6f}"+'，峰對峰比為 '+f"{ratios[0]['ptp_ratio']:.6f}"+'。這是完整臂別記錄中混合活動比例的比較，不是單次最大收縮、肌肉力量或相同負載下的效果比。兩段時長不同（18.468、17.868 秒），RMS 已按樣本數正規化，但收縮占比仍可能不同，峰對峰也受記錄長度與孤立尖波影響。')
    fig('EMG_full_channels_events','圖 6｜實際 EMG 原始與積分通道。18.468 秒為主力臂／非主力臂切換標記；沒有單次握緊起訖標記。')
    fig('EMG_arm_comparison','圖 7｜完整標記臂別區段的去平均 RMS 與峰對峰；含活動與低活動間隔。')
    p('未完成的條件比較：放鬆對用力、逐次握力或不同負載，因為沒有可重現的實際狀態起訖標記。內嵌 Journal 雖列四次握緊流程與既有 mV-sec 表值，並未給出每次量測的時間窗；既有表值不等於本次 CH1 去平均 RMS，未拿來代填。完整標記臂別比較已完成，沒有自行把第四個活動群命名為最大握力。')
    h('結果適用範圍與限制')
    p('ECG：有實際條件標記且 R 峰可辨識，因此可完成心率與 R–R 描述性比較。基線位移不影響所有 R 波的可見性，但限制精細形態判讀；此次沒有計算 PR、QRS、QT、ST 或診斷性指標。')
    p('EMG：可辨識臂別並計算完整區段 RMS／峰對峰；不能把未標記的狀態當成已確認的用力、放鬆或負載條件。原始積分通道保持 mV-sec，不拿來替代本次 RMS。')
    p('這些是單次記錄中的描述。相鄰 R–R、呼吸回合或 EMG 樣本不當作獨立受試者，沒有群體顯著性檢定。使用者已確認這是 HW1 的直接 BIOPAC 記錄，作為 HW2 輸入。先前 EEG 可能有一份作廢的提醒不自動套用到此次 ECG／EMG，這兩份也沒有依預期結果篩選或宣告作廢。')
    h('可重現程式與交付')
    p('執行：python -m pip install -r requirements.txt；python analyze_ecg_emg.py --input raw --output results。程式依實際通道名辨識 ECG／EMG，並檢查條件標記、校正、取樣率及 R–R 一致性；原始輸入 SHA-256 符合本次資料時，會重建這份含人工目視註解的中文報告。')
    p('版本：Python '+platform.python_version()+'；'+'；'.join(k+' '+v for k,v in versions.items())+'。主程式、精確版本需求、未修改原始資料與全部 CSV／PNG／SVG 均包含在完整包中。')
    table(['檔案','內容'],[
        ['ECG_condition_summary.csv','各實際條件的區段、R點數、R–R及心率'],
        ['ECG_R_peaks.csv / ECG_RR_intervals.csv','每個偵測 R 點及逐間期結果，可核對與重畫'],
        ['ECG_detector_sensitivity.csv','0.10、0.12、0.15 mV 的候選點比較'],
        ['ECG_breath_phase_summary.csv / ECG_breath_phase_RR.csv','依實際呼吸標記界定的補充摘要'],
        ['EMG_arm_summary.csv / EMG_arm_ratios.csv','只含可辨識臂別條件的原始EMG統計'],
        ['ECG_full_channels.csv / EMG_full_channels.csv','完整已校正通道序列；CH40保留其原始單位'],
        ['channel_metadata.csv / event_markers.csv','實際名稱、單位、fs、SHA-256與全部原始標記'],
        ['quality_summary.csv / quality_1s_blocks.csv / quality_visual_annotations.csv','數值品質指標、1秒統計與目視註解'],
        ['analysis_manifest.json / validation_checks.json','可重現設定、版本及數值一致性核對'],
    ])
    fig('ECG_R_peak_full_review','補圖｜全程 R 波疊圖，每面20秒、各面獨立縱軸；圖中包含分析首尾1秒之外的偵測點以供核對。')
    css='body{font-family:system-ui,"Microsoft JhengHei","Noto Sans CJK TC",sans-serif;max-width:1200px;margin:36px auto;padding:0 24px;color:#203247;line-height:1.8}h1{font-size:29px;line-height:1.4}h2{margin-top:36px;font-size:22px;color:#165477;border-bottom:1px solid #c8dbe4;padding-bottom:8px}p{margin:14px 0}.table-wrap{overflow:auto;margin:18px 0}table{width:100%;border-collapse:collapse;font-size:14px;line-height:1.6}th,td{padding:10px;border:1px solid #d7e1e8;text-align:left;vertical-align:top}th{background:#eaf1f6}tr:nth-child(even){background:#f8fafc}figure{margin:26px 0}img{width:100%;height:auto;border:1px solid #dce4ea}figcaption{font-size:14px;color:#526574}@media print{body{max-width:none;margin:0;font-size:10pt}h2{break-after:avoid}figure,tr{break-inside:avoid}.table-wrap{overflow:visible}}'
    (out.parent/'Group2_ECG_EMG_Report.html').write_text('<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Group 2 ECG／EMG 分析</title><style>'+css+'</style><body>'+''.join(blocks)+'</body></html>',encoding='utf-8')
    (out.parent/'Group2_ECG_EMG_結果與限制.md').write_text('\n'.join(md),encoding='utf-8')


def run(root,out):
    out.mkdir(parents=True,exist_ok=True)
    data={};metadata=[];events=[];quality=[];blocks=[];checks=[]
    candidates=[]
    for p in root.rglob('*'):
        if p.is_file() and (p.name.endswith('-L05') or p.name.endswith('-L01')):candidates.append(p)
    for p in sorted(candidates):
        d=bioread.read_file(str(p));ch=d.channels[0]
        if ch.name not in ['ECG','EMG']:raise ValueError('Unexpected first-channel name')
        label=ch.name
        if label in data:raise ValueError('Duplicate input modality')
        data[label]=d;fs=ch.samples_per_second
        if ch.units!='mV':raise ValueError('Unexpected main channel unit; no implicit conversion')
        assert np.allclose(ch.data,ch.raw_data.astype(float)*ch.raw_scale_factor+ch.raw_offset)
        digest=hashlib.sha256(p.read_bytes()).hexdigest()
        for c in d.channels:
            metadata.append(dict(file=label,filename=p.name,sha256=digest,channel=c.order_num,name=c.name,
                                 original_unit=c.units,fs_Hz=c.samples_per_second,samples=c.point_count,
                                 duration_s=c.point_count/c.samples_per_second,
                                 last_sample_time_s=(c.point_count-1)/c.samples_per_second,
                                 storage_dtype=str(c.dtype),raw_scale=c.raw_scale_factor,raw_offset=c.raw_offset,
                                 file_source='BIOPAC/BSL native record',course_exp_number='HW1, reused as HW2 input (user confirmed)'))
        for e in d.event_markers:
            events.append(dict(file=label,sample_index_zero_based=e.sample_index,time_s=e.time_index,
                               original_text=e.text,type_code=e.type_code,
                               marker_role='automatic rate output' if e.type_code=='nrto' else 'recorded condition/phase'))
        quality.append(qc(label,'Full record',ch,0,len(ch.data)))
        for i in range(0,len(ch.data),round(fs)):
            v=ch.data[i:i+round(fs)]
            blocks.append(dict(file=label,start_s=i/fs,end_exclusive_s=(i+len(v))/fs,samples=len(v),
                               mean_mV=float(v.mean()),demeaned_rms_mV=float(v.std()),ptp_mV=float(np.ptp(v))))
        headers=['time_s']+[f'ch{c.order_num}_{"mV" if c.order_num==1 else "native_value"}' for c in d.channels]
        assert all(c.samples_per_second==fs and len(c.data)==len(ch.data) for c in d.channels)
        write_csv(out/f'{label}_full_channels.csv',
                  (dict(zip(headers,[i/fs]+[float(c.data[i]) for c in d.channels])) for i in range(len(ch.data))))
    if set(data)!= {'ECG','EMG'}:raise ValueError('Need actual ECG and EMG inputs')

    # ECG: use only recorded apnd conditions; no intervals cross append boundaries.
    d=data['ECG'];ch=d.channels[0];x=ch.data;fs=ch.samples_per_second;t=np.arange(len(x))/fs
    conditions=[e for e in d.event_markers if e.type_code=='apnd']
    if [e.text.split(',')[0] for e in conditions] != ECG_ZH:raise ValueError('Different ECG condition markers')
    if fs!=1000:raise ValueError('Unexpected ECG fs')
    peaks_rows=[];rr_rows=[];ecg_summary=[];sensitivity=[];condition_ranges=[];all_peaks=[];validation_peaks=[]
    for i,e in enumerate(conditions):
        a=e.sample_index;b=conditions[i+1].sample_index if i+1<len(conditions) else len(x)
        guard=round(DETECTOR['analysis_boundary_guard_s']*fs);aa=a+guard;bb=b-guard
        condition_ranges.append((aa,bb))
        quality.append(qc('ECG',ECG_ZH[i]+' / analysis',ch,aa,bb))
        baseline,_=r_peaks(x[a:b],fs,.12)
        baseline=baseline+a;all_peaks.extend(baseline)
        accepted=baseline[(baseline>=aa)&(baseline<bb)]
        for th in [.10,.12,.15]:
            pp,_=r_peaks(x[a:b],fs,th);pp=pp+a;pp=pp[(pp>=aa)&(pp<bb)]
            sensitivity.append(dict(condition=ECG_ZH[i],prominence_mV=th,peaks=len(pp),
                                    identical_times_to_0p12=bool(np.array_equal(accepted,pp)),
                                    extra_peaks=int(len(np.setdiff1d(pp,accepted))),
                                    missing_peaks=int(len(np.setdiff1d(accepted,pp)))))
            if th==.15:validation_peaks.append(bool(np.array_equal(accepted,pp)))
        for j,p in enumerate(baseline):
            peaks_rows.append(dict(condition=ECG_ZH[i],r_index_zero_based=int(p),r_time_s=p/fs,
                                   raw_R_amplitude_mV=float(x[p]),
                                   included_in_summary=bool(aa<=p<bb),
                                   inclusion_reason='inside guarded condition' if aa<=p<bb else '1s boundary guard',
                                   detector_prominence_mV=.12))
        rr=np.diff(accepted)/fs
        if not len(rr):raise ValueError('Insufficient R peaks')
        for j,(p,q,interval) in enumerate(zip(accepted[:-1],accepted[1:],rr)):
            rr_rows.append(dict(condition=ECG_ZH[i],rr_number_in_condition=j+1,
                                r1_index=int(p),r2_index=int(q),r1_time_s=p/fs,r2_time_s=q/fs,
                                midpoint_s=(p+q)/(2*fs),RR_s=float(interval),RR_ms=float(interval*1000),
                                instantaneous_HR_bpm=float(60/interval),
                                rhythm_review_flag=bool(interval<.5*np.median(rr) or interval>1.5*np.median(rr))))
        ecg_summary.append(dict(condition=ECG_ZH[i],record_start_s=a/fs,record_end_exclusive_s=b/fs,
                                analysis_start_index=aa,analysis_stop_index_exclusive=bb,
                                analysis_start_s=aa/fs,analysis_end_exclusive_s=bb/fs,
                                analysis_duration_s=(bb-aa)/fs,R_peaks=len(accepted),RR_intervals=len(rr),
                                first_R_s=accepted[0]/fs,last_R_s=accepted[-1]/fs,
                                RR_mean_ms=float(rr.mean()*1000),RR_std_ms_ddof1=float(rr.std(ddof=1)*1000),
                                RR_min_ms=float(rr.min()*1000),RR_max_ms=float(rr.max()*1000),
                                HR_from_mean_RR_bpm=float(60/rr.mean()),
                                mean_instantaneous_HR_bpm=float(np.mean(60/rr)),
                                HR_min_bpm=float(np.min(60/rr)),HR_max_bpm=float(np.max(60/rr)),
                                native_CH40_mean_bpm=float(d.channels[1].data[aa:bb].mean()),
                                waveform_scope='R timing accepted after visual overlay review; baseline drift limits morphology analysis'))
        checks.append(dict(file='ECG',condition=ECG_ZH[i],calibration_check=True,
                           RR_count_equals_Rminus1=(len(rr)==len(accepted)-1),
                           no_RR_crosses_condition_boundary=bool(np.all(accepted>=aa)&np.all(accepted<bb)),
                           sum_RR_equals_last_minus_first=bool(np.isclose(rr.sum(),(accepted[-1]-accepted[0])/fs)),
                           manual_R_additions_or_deletions=0))

    all_peaks=np.array(all_peaks)
    write_csv(out/'ECG_R_peaks.csv',peaks_rows);write_csv(out/'ECG_RR_intervals.csv',rr_rows)
    write_csv(out/'ECG_condition_summary.csv',ecg_summary)
    write_csv(out/'ECG_detector_sensitivity.csv',sensitivity)
    assert all(validation_peaks),'0.12 vs 0.15 mV mismatch: inspect before accepting'

    # Deep-breath phase summaries: recorded markers only, no RR straddles a phase.
    phase_markers=[e for e in d.event_markers if e.type_code=='defl']
    breath_rows=[];phase_intervals=[]
    for j,e in enumerate(phase_markers):
        a=max(e.sample_index,condition_ranges[2][0])
        b=min(phase_markers[j+1].sample_index if j+1<len(phase_markers) else conditions[3].sample_index,condition_ranges[2][1])
        phase=e.text
        applicable=[r for r in rr_rows if r['condition']=='深呼吸' and r['r1_index']>=a and r['r2_index']<b]
        for r in applicable:phase_intervals.append(dict(phase=phase,phase_number=j+1,**r))
        breath_rows.append(dict(phase_number=j+1,recorded_phase_marker=phase,marker_time_s=e.time_index,
                               analysis_phase_start_s=a/fs,analysis_phase_end_exclusive_s=b/fs,
                               RR_intervals=len(applicable),
                               RR_mean_ms=float(np.mean([r['RR_ms'] for r in applicable])) if applicable else '',
                               HR_from_mean_RR_bpm=60/np.mean([r['RR_s'] for r in applicable]) if applicable else '',
                               scope='Whole RR intervals contained in this recorded phase only; descriptive'))
    write_csv(out/'ECG_breath_phase_summary.csv',breath_rows)
    write_csv(out/'ECG_breath_phase_RR.csv',phase_intervals)

    fig,axs=plt.subplots(3,1,figsize=(14,10),sharex=True,layout='constrained')
    axs[0].plot(t,x,lw=.5,color='#253c54');axs[0].set_ylabel('Raw ECG (mV)')
    axs[1].plot(t,d.channels[1].data,color='#999999',lw=.7,label='Native CH40 (not accepted as primary HR)')
    for i,condition in enumerate(ECG_ZH):
        rr=[r for r in rr_rows if r['condition']==condition]
        axs[1].plot([r['midpoint_s'] for r in rr],[r['instantaneous_HR_bpm'] for r in rr],'.-',ms=3,lw=.8,color=COLORS[i],label=ECG_LABELS[i]+' / redetected')
    axs[1].set_ylabel('Heart rate (bpm)');axs[1].legend(fontsize=8,ncol=2)
    br=[b for b in blocks if b['file']=='ECG']
    axs[2].plot([b['start_s'] for b in br],[b['mean_mV'] for b in br],lw=1,color='#647b91')
    axs[2].set_ylabel('Nonoverlap 1-s mean (mV)')
    for ax in axs:decorate_events(ax,conditions,fs,len(x)/fs,ECG_LABELS,annotate=(ax is axs[0]))
    axs[-1].set_xlabel('File time (s); appended records, not necessarily continuous wall-clock time')
    fig.suptitle('ECG / recorded conditions / full waveform and native-vs-redetected HR')
    savefig(fig,out,'ECG_full_events_quality')

    # All samples and all detected R peaks can be visually reviewed in 20-s panels.
    fig,axs=plt.subplots(8,1,figsize=(14,18),layout='constrained')
    for j,ax in enumerate(axs):
        a=j*20000;b=min((j+1)*20000,len(x));pp=all_peaks[(all_peaks>=a)&(all_peaks<b)]
        ax.plot(t[a:b],x[a:b],lw=.55,color='#436982');ax.plot(pp/fs,x[pp],'o',color='#c7354d',ms=2.6)
        ax.set(title=f'{a/fs:.0f}–{b/fs:.3f} s',ylabel='mV',xlim=(a/fs,b/fs));ax.grid(alpha=.18)
        for e in conditions:
            if a/fs<=e.time_index<b/fs:ax.axvline(e.time_index,color='#333333',ls='--',lw=.8)
    axs[-1].set_xlabel('Time (s); red dots = detected R; panels have independent y axes')
    savefig(fig,out,'ECG_R_peak_full_review')

    fig,axs=plt.subplots(4,1,figsize=(13,10),layout='constrained')
    for i,(ax,(a,b)) in enumerate(zip(axs,[(2,10),(23,31),(45,53),(87,95)])):
        aa=round(a*fs);bb=round(b*fs);pp=all_peaks[(all_peaks>=aa)&(all_peaks<bb)]
        ax.plot(t[aa:bb],x[aa:bb],lw=.75,color=COLORS[i]);ax.plot(pp/fs,x[pp],'o',ms=4,color='#c7354d')
        ax.set(title=ECG_LABELS[i]+' / R-wave timing review',ylabel='mV',xlim=(a,b));ax.grid(alpha=.18)
    axs[-1].set_xlabel('File time (s)');savefig(fig,out,'ECG_R_peak_examples')

    fig,axs=plt.subplots(2,1,figsize=(14,7),sharex=True,layout='constrained')
    for i,condition in enumerate(ECG_ZH):
        rr=[r for r in rr_rows if r['condition']==condition]
        axs[0].plot([r['midpoint_s'] for r in rr],[r['RR_ms'] for r in rr],'.-',color=COLORS[i],ms=3,lw=.8,label=ECG_LABELS[i])
        axs[1].plot([r['midpoint_s'] for r in rr],[r['instantaneous_HR_bpm'] for r in rr],'.-',color=COLORS[i],ms=3,lw=.8)
    for ax in axs:
        decorate_events(ax,conditions,fs,len(x)/fs,ECG_LABELS,annotate=False)
        for e in phase_markers:ax.axvline(e.time_index,color='#93a7af',lw=.5,alpha=.7)
    axs[0].set_ylabel('R–R interval (ms)');axs[1].set(ylabel='60 / R–R (bpm)',xlabel='File time (s)')
    axs[0].legend(ncol=4,fontsize=9)
    fig.suptitle('ECG / only intervals fully inside each guarded condition; no RR deletion')
    savefig(fig,out,'ECG_RR_HR_timeline')

    fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
    for ax,key,title,unit in [(axs[0],'HR_from_mean_RR_bpm','Heart rate = 60 / mean R–R','bpm'),
                             (axs[1],'RR_mean_ms','Mean R–R interval','ms')]:
        vals=[r[key] for r in ecg_summary];ax.bar(ECG_LABELS,vals,color=COLORS)
        for j,v in enumerate(vals):ax.text(j,v,f'{v:.2f}',ha='center',va='bottom')
        ax.set(title=title,ylabel=unit,ylim=(0,max(vals)*1.18));ax.tick_params(axis='x',rotation=15);ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
    fig.suptitle('ECG / one record per condition; descriptive values, no inferential error bars')
    savefig(fig,out,'ECG_condition_comparison')

    # Supine detector discrepancy; native rate-update marker times are retained.
    fig,axs=plt.subplots(2,1,figsize=(12,6),layout='constrained')
    a,b=2000,6000;pp=all_peaks[(all_peaks>=a)&(all_peaks<b)]
    axs[0].plot(t[a:b],x[a:b],lw=.8,label='Raw ECG');axs[0].plot(pp/fs,x[pp],'o',ms=5,color='#c7354d',label='Redetected R')
    nrto=[e for e in d.event_markers if e.type_code=='nrto' and a/fs<=e.time_index<b/fs]
    for e in nrto:axs[0].axvline(e.time_index,color='#daa045',alpha=.6,lw=.9)
    axs[0].set(xlim=(2,6),ylabel='mV',title='Supine: native HR output disagrees with raw R–R timing');axs[0].legend()
    axs[1].plot(t[a:b],d.channels[1].data[a:b],lw=1,color='#777777',label='Native CH40')
    rr=[r for r in rr_rows if r['condition']=='仰臥' and 2<=r['midpoint_s']<6]
    axs[1].plot([r['midpoint_s'] for r in rr],[r['instantaneous_HR_bpm'] for r in rr],'o-',color='#c7354d',label='Redetected 60/RR')
    axs[1].set(xlim=(2,6),ylabel='bpm',xlabel='File time (s)');axs[1].legend()
    for ax in axs:ax.grid(alpha=.18)
    savefig(fig,out,'ECG_native_HR_discrepancy')

    # EMG: only the two actual arm markers define conditions. No force/load labels.
    d=data['EMG'];ch=d.channels[0];x=ch.data;fs=ch.samples_per_second;t=np.arange(len(x))/fs
    arms=[e for e in d.event_markers if e.type_code=='apnd']
    if [e.text.split(',')[0] for e in arms]!=EMG_ZH:raise ValueError('Different EMG arm markers')
    if fs!=500:raise ValueError('Unexpected EMG fs')
    emg_summary=[]
    for i,e in enumerate(arms):
        a=e.sample_index;b=arms[i+1].sample_index if i+1<len(arms) else len(x);v=x[a:b]
        q=qc('EMG',EMG_ZH[i],ch,a,b);quality.append(q)
        emg_summary.append(dict(condition=EMG_ZH[i],start_index=a,stop_index_exclusive=b,
                                start_s=a/fs,end_exclusive_s=b/fs,samples=len(v),duration_s=len(v)/fs,
                                mean_mV=float(v.mean()),demeaned_rms_mV=float(np.sqrt(np.mean((v-v.mean())**2))),
                                ptp_mV=float(np.ptp(v)),min_mV=float(v.min()),max_mV=float(v.max()),
                                condition_basis='actual dominant/nondominant arm event marker',
                                contraction_relaxation_or_load='not separately marked; mixed activity period',
                                scope='Whole marked arm interval; not a matched contraction or maximal-force comparison'))
        checks.append(dict(file='EMG',condition=EMG_ZH[i],calibration_check=True,
                           RMS_equals_population_std=bool(np.isclose(np.sqrt(np.mean((v-v.mean())**2)),v.std(ddof=0))),
                           segment_length_matches_samples=(len(v)==b-a),no_samples_removed=True))
    ratios=[dict(comparison='非主力臂 / 主力臂（完整標記區段）',
                 demeaned_rms_ratio=emg_summary[1]['demeaned_rms_mV']/emg_summary[0]['demeaned_rms_mV'],
                 ptp_ratio=emg_summary[1]['ptp_mV']/emg_summary[0]['ptp_mV'],
                 scope='Different full-record durations and mixed duty cycles; not muscle force ratio')]
    write_csv(out/'EMG_arm_summary.csv',emg_summary);write_csv(out/'EMG_arm_ratios.csv',ratios)
    fig,axs=plt.subplots(2,1,figsize=(14,7),sharex=True,layout='constrained')
    axs[0].plot(t,x,lw=.5,color='#245d89');axs[0].set_ylabel('Raw EMG (mV)')
    axs[1].plot(t,d.channels[1].data,lw=.7,color='#556e80');axs[1].set_ylabel('Native integrated EMG (mV-sec)')
    for ax in axs:decorate_events(ax,arms,fs,len(x)/fs,EMG_LABELS,annotate=(ax is axs[0]))
    axs[1].set_xlabel('File time (s); individual contractions and loads are not event-marked')
    fig.suptitle('EMG (actual channel) / filename contains ecg / recorded arm conditions only')
    savefig(fig,out,'EMG_full_channels_events')
    fig,axs=plt.subplots(1,2,figsize=(10,4.6),layout='constrained')
    for ax,key,title in [(axs[0],'demeaned_rms_mV','Demeaned RMS'),(axs[1],'ptp_mV','Peak-to-peak')]:
        vals=[r[key] for r in emg_summary];ax.bar(EMG_LABELS,vals,color=COLORS[:2],width=.6)
        for j,v in enumerate(vals):ax.text(j,v,f'{v:.6f}',ha='center',va='bottom')
        ax.set(title=title,ylabel='mV',ylim=(0,max(vals)*1.2));ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
    fig.suptitle('EMG / whole marked arm intervals (18.468 vs 17.868 s; mixed activity)')
    savefig(fig,out,'EMG_arm_comparison')

    write_csv(out/'channel_metadata.csv',metadata);write_csv(out/'event_markers.csv',events)
    write_csv(out/'quality_summary.csv',quality);write_csv(out/'quality_1s_blocks.csv',blocks)
    versions={k:importlib.metadata.version(k) for k in ['numpy','scipy','matplotlib','bioread']}
    manifest=dict(python=platform.python_version(),versions=versions,detector=DETECTOR,
                  ECG_boundary_rule='1s guard at each appended condition start/end; no cross-boundary RR',
                  EMG_boundary_rule='Full actual arm event intervals; no contraction state inference',
                  signal_scaling='bioread calibrated data once; no additional gain scaling',
                  source='BIOPAC/BSL HW1 recordings; user confirmed direct electrode-to-BIOPAC acquisition and reuse as HW2 inputs; no additional participant acquisition required',
                  no_manual_R_edits=True,no_RR_deletion=True,no_group_significance_testing=True,
                  discarded_record_identity='Not established for these files; no prior EEG warning transferred',
                  channels=metadata,validation_checks=checks,stable_0p12_vs_0p15_mV=validation_peaks)
    (out/'analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'validation_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
    if all(next(m['sha256'] for m in metadata if m['file']==label)==h for label,h in KNOWN.items()):
        make_report(out,metadata,events,quality,ecg_summary,emg_summary,ratios,breath_rows,sensitivity,versions)
    else:print('Different input hashes: dataset-specific manual QC report not applied.')
    print(json.dumps(dict(ECG=ecg_summary,EMG=emg_summary,EMG_ratios=ratios),ensure_ascii=False,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,default=Path(__file__).parent/'raw')
    p.add_argument('--output',type=Path,default=Path(__file__).parent/'results')
    args=p.parse_args();run(args.input,args.output)


if __name__=='__main__':
    main()
