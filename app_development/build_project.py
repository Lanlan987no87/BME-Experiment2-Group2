"""Modify the user's working AI2 project; keep BLE extension binaries unchanged."""
from pathlib import Path
import json,zipfile,itertools,xml.etree.ElementTree as E,hashlib
R=Path(__file__).resolve().parent;BASE=R.parent;SRC=BASE/'upload/Group2_Biosignal.aia';OUT=R/'Group2_Biosignal.aia'
ids=itertools.count(1)
def b(kind,**atts):return E.Element('block',dict(type=kind,id='g2_'+str(next(ids)),**{k:str(v) for k,v in atts.items()}))
def field(n,name,value):E.SubElement(n,'field',name=name).text=str(value);return n
def val(n,name,value):E.SubElement(n,'value',name=name).append(value);return n
def mut(n,**atts):return E.SubElement(n,'mutation',dict(xmlns='http://www.w3.org/1999/xhtml',**{k:str(v) for k,v in atts.items()}))
def chain(*blocks):
 blocks=[x for x in blocks if x is not None]
 for a,z in zip(blocks,blocks[1:]):E.SubElement(a,'next').append(z)
 return blocks[0] if blocks else None
def stmt(n,name,*blocks):
 c=chain(*blocks)
 if c is not None:E.SubElement(n,'statement',name=name).append(c)
 return n
def text(s):return field(b('text'),'TEXT',s)
def num(x):return field(b('math_number'),'NUM',x)
def boolean(x):return field(b('logic_boolean'),'BOOL','TRUE' if x else 'FALSE')
def get(name,event=False):
 n=b('lexical_variable_get')
 if event:E.SubElement(mut(n),'eventparam',name=name)
 return field(n,'VAR',name)
def gg(name):return get('global '+name)
def gs(name,v):return val(field(b('lexical_variable_set'),'VAR','global '+name),'VALUE',v)
def prop(inst,typ,name,v=None):
 n=b('component_set_get');mut(n,component_type=typ,set_or_get='get' if v is None else 'set',property_name=name,is_generic='false',instance_name=inst);field(n,'COMPONENT_SELECTOR',inst);field(n,'PROP',name)
 return val(n,'VALUE',v) if v is not None else n
def method(inst,typ,name,*args):
 n=b('component_method');mut(n,component_type=typ,method_name=name,is_generic='false',instance_name=inst);field(n,'COMPONENT_SELECTOR',inst)
 for i,a in enumerate(args):val(n,'ARG'+str(i),a)
 return n
def event(inst,typ,name,*body,x=0,y=0):
 n=b('component_event',x=x,y=y);mut(n,component_type=typ,is_generic='false',instance_name=inst,event_name=name);field(n,'COMPONENT_SELECTOR',inst);return stmt(n,'DO',*body)
def join(*xs):
 n=b('text_join');mut(n,items=len(xs))
 for i,x in enumerate(xs):val(n,'ADD'+str(i),x)
 return n
def eq(a,z):return val(val(field(b('logic_compare'),'OP','EQ'),'A',a),'B',z)
def cmp(a,z,op):return val(val(field(b('math_compare'),'OP',op),'A',a),'B',z)
def add(a,z):
 n=b('math_add');mut(n,items=2);return val(val(n,'NUM0',a),'NUM1',z)
def length(x):return val(b('text_length'),'VALUE',x)
def iff(cond,then,otherwise=None):
 n=b('controls_if')
 if otherwise:mut(n,**{'else':'1'})
 val(n,'IF0',cond);stmt(n,'DO0',*then)
 if otherwise:stmt(n,'ELSE',*otherwise)
 return n
def status(s):return [gs('Status',text(s)),prop('LblStatus','Label','Text',text(s))]
def buttons(scan,stop,disc):return [prop(i,'Button','Enabled',boolean(v)) for i,v in [('BtnScan',scan),('BtnStopScan',stop),('BtnDisconnect',disc)]]
def stopped():return [method('BluetoothLE1','BluetoothLE','StopScanning'),prop('ScanClock','Clock','TimerEnabled',boolean(False))]
def base_state(s):return [*stopped(),*status(s),*buttons(True,False,False),prop('ListDevices','ListView','Visible',boolean(True))]
root=E.Element('xml',xmlns='https://developers.google.com/blockly/xml')
for i,(name,v) in enumerate([('Rx',text('')),('Session',num(0)),('Overflow',num(0)),('Status',text('尚未連線'))]):
 n=b('global_declaration',x=20,y=20+i*80);field(n,'NAME',name);val(n,'VALUE',v);root.append(n)
root.append(event('Screen1','Form','Initialize',prop('WebViewer1','WebViewer','WebViewString',text('')),*status('尚未連線'),*buttons(True,False,False),x=20,y=390))
root.append(event('BtnScan','Button','Click',*status('掃描中，請選 HMSoft'),prop('ListDevices','ListView','Visible',boolean(True)),prop('ListDevices','ListView','SelectionIndex',num(0)),*buttons(False,True,False),method('BluetoothLE1','BluetoothLE','StartScanning'),prop('ScanClock','Clock','TimerEnabled',boolean(True)),x=500,y=20))
root.append(event('BtnStopScan','Button','Click',*stopped(),*status('掃描已停止；可選裝置或重新掃描'),*buttons(True,False,False),x=500,y=360))
root.append(event('ScanClock','Clock','Timer',*stopped(),*status('掃描完成，請選 HMSoft'),*buttons(True,False,False),x=500,y=610))
root.append(event('BluetoothLE1','BluetoothLE','DeviceFound',prop('ListDevices','ListView','ElementsFromString',prop('BluetoothLE1','BluetoothLE','DeviceList')),x=500,y=850))
root.append(event('ListDevices','ListView','AfterPicking',iff(cmp(prop('ListDevices','ListView','SelectionIndex'),num(0),'GT'),[*stopped(),*status('連線中…'),*buttons(False,False,False),prop('ListDevices','ListView','Visible',boolean(False)),method('BluetoothLE1','BluetoothLE','Connect',prop('ListDevices','ListView','SelectionIndex'))]),x=1050,y=20))
root.append(event('BluetoothLE1','BluetoothLE','Connected',gs('Rx',text('')),gs('Session',add(gg('Session'),num(1))),*status('已連線'),*buttons(False,False,True),prop('ListDevices','ListView','Visible',boolean(False)),method('BluetoothLE1','BluetoothLE','RegisterForStrings',text('0000ffe0-0000-1000-8000-00805f9b34fb'),text('0000ffe1-0000-1000-8000-00805f9b34fb'),boolean(False)),x=1050,y=440))
root.append(event('BtnDisconnect','Button','Click',method('BluetoothLE1','BluetoothLE','Disconnect'),x=1050,y=940))
root.append(event('BluetoothLE1','BluetoothLE','Disconnected',gs('Rx',text('')),*base_state('已中斷連線'),x=1600,y=20))
root.append(event('BluetoothLE1','BluetoothLE','ConnectionFailed',*base_state('連線失敗，請重新掃描'),method('Notifier1','Notifier','ShowAlert',join(text('連線失敗：'),get('reason',True))),x=1600,y=380))
loop=b('controls_forEach');field(loop,'VAR','item');val(loop,'LIST',get('stringValues',True));stmt(loop,'DO',gs('Rx',join(gg('Rx'),get('item'))))
# Explicit bounded native queue. Overflow discards that batch and resets the JS partial frame through the counter.
root.append(event('BluetoothLE1','BluetoothLE','StringsReceived',iff(prop('BluetoothLE1','BluetoothLE','IsDeviceConnected'),[loop,iff(cmp(length(gg('Rx')),num(32768),'GT'),[gs('Rx',text('')),gs('Overflow',add(gg('Overflow'),num(1)))])]),x=2200,y=20))
# JavaScript requests one batch; native replies once. No overwrite of an unconsumed envelope.
starts=val(val(b('text_starts_at'),'TEXT',get('value',True)),'PIECE',text('CSV:'))
sub=val(val(val(b('text_segment'),'TEXT',get('value',True)),'START',num(5)),'LENGTH',val(val(b('math_subtract'),'A',length(get('value',True))),'B',num(4)))
csvif=iff(eq(starts,num(1)),[method('File1','File','SaveFile',sub,text('Group2_ADC_latest100.csv')),prop('WebViewer1','WebViewer','WebViewString',text('ACK'))])
root.append(event('WebViewer1','WebViewer','WebViewStringChange',iff(eq(get('value',True),text('POLL')),[prop('WebViewer1','WebViewer','WebViewString',join(gg('Session'),text('|'),gg('Overflow'),text('|'),gg('Status'),text('|'),gg('Rx'))),gs('Rx',text(''))],[csvif]),x=2200,y=650))
root.append(event('File1','File','AfterFileSaved',method('Sharing1','Sharing','ShareFile',get('fileName',True)),x=3000,y=20))
root.append(event('Screen1','Form','ErrorOccurred',method('Notifier1','Notifier','ShowAlert',join(text('操作訊息：'),get('message',True))),x=3000,y=200))
root.append(event('Screen1','Form','PermissionDenied',method('Notifier1','Notifier','ShowAlert',join(text('請在 Android 設定允許附近裝置／藍牙權限：'),get('permissionName',True))),x=3000,y=410))
E.SubElement(root,'yacodeblocks',{'xmlns':'https://appinventor.mit.edu/ns/project/','ya-version':'237','language-version':'39'})
# Designer components; original BLE extension and account namespace remain intact.
uids=itertools.count(2000000000)
def comp(name,typ,ver,**props):return dict({'$Name':name,'$Type':typ,'$Version':str(ver),'Uuid':str(next(uids))},**{k:str(v) for k,v in props.items()})
with zipfile.ZipFile(SRC) as z:
 entries={n:z.read(n) for n in z.namelist()}
scmname=next(n for n in entries if n.endswith('Screen1.scm'));bkyname=scmname[:-3]+'bky'
scm=json.loads(entries[scmname].decode().split('$JSON\n',1)[1].rsplit('\n|#',1)[0]);screen=scm['Properties'];oldble=next(c for c in screen['$Components'] if c['$Type']=='BluetoothLE')
oldble.update(AutoReconnect='False',ConnectionTimeout='10')
main=comp('MainLayout','VerticalArrangement',4,Width=-2,Height=-2)
bar=comp('ConnectionBar','HorizontalArrangement',4,Width=-2)
bar['$Components']=[comp('BtnScan','Button',8,Text='掃描藍牙',FontSize=14),comp('BtnStopScan','Button',8,Text='停止掃描',Enabled='False',FontSize=14),comp('BtnDisconnect','Button',8,Text='中斷連線',Enabled='False',FontSize=14)]
main['$Components']=[comp('LblTitle','Label',6,Text='Group 2 生醫訊號',FontSize=22,FontBold='True',Width=-2),bar,comp('LblStatus','Label',6,Text='尚未連線',FontSize=15,Width=-2),comp('ListDevices','ListView',11,Width=-2,Height=120,FontSize=15,BackgroundColor='&HFFFFFFFF',TextColor='&HFF17394D'),comp('WebViewer1','WebViewer',11,Width=-2,Height=-2,HomeUrl='http://localhost/waveform.html',UsesLocation='False',FollowLinks='False')]
screen.update(Title='Group2 Biosignal',AppName='Group2 Biosignal',Scrollable='False',Sizing='Responsive',ScreenOrientation='unspecified',BackgroundColor='&HFFF1F5F9')
screen['$Components']=[main,oldble,comp('ScanClock','Clock',4,TimerEnabled='False',TimerInterval=10000,TimerAlwaysFires='False'),comp('File1','File',4,DefaultScope='App'),comp('Sharing1','Sharing',2),comp('Notifier1','Notifier',7)]
entries[scmname]=('#|\n$JSON\n'+json.dumps(scm,ensure_ascii=True,separators=(',',':'))+'\n|#\n').encode()
entries[bkyname]=E.tostring(root,encoding='utf-8')
entries['assets/waveform.html']=(R/'waveform.html').read_bytes()
propname='youngandroidproject/project.properties';s=entries[propname].decode().replace('versioncode=1','versioncode=2').replace('versionname=1.0','versionname=2.0').replace('androidminsdk=14','androidminsdk=21');entries[propname]=s.encode()
with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED) as z:
 for n,data in entries.items():z.writestr(n,data)
print(OUT,OUT.stat().st_size,'bytes;',sum(1 for _ in root.iter('block')),'blocks')
# Structural integrity and exactly-preserved extension byte checks.
component_names={}
def visit(c):
 component_names[c['$Name']]=c['$Type']
 for child in c.get('$Components',[]):visit(child)
visit(screen)
for m in root.iter('mutation'):
 if 'instance_name' in m.attrib:assert component_names[m.attrib['instance_name']]==m.attrib['component_type']
allids=[n.attrib['id'] for n in root.iter('block')];assert len(allids)==len(set(allids))
with zipfile.ZipFile(OUT) as z:
 assert z.testzip() is None
 for n in entries:
  if n.startswith('assets/external_comps/'):assert z.read(n)==zipfile.ZipFile(SRC).read(n)
(R/'project_validation.json').write_text(json.dumps({'block_count':len(allids),'component_count':len(component_names),'extension_binaries_unchanged':True,'zip_crc_check':'passed','xml_json_parse':'passed','hardware_test':'pending on user tablet','apk_build':'not performed; use MIT build menu'},indent=2))
