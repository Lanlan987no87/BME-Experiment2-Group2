"""Execute generated Blockly control flow against mock components; no Android claim."""
from pathlib import Path
import zipfile,json,xml.etree.ElementTree as E
R=Path(__file__).resolve().parent;z=zipfile.ZipFile(R/'Group2_Biosignal.aia')
r=E.fromstring(z.read(next(n for n in z.namelist() if n.endswith('.bky'))))
for el in r.iter():el.tag=el.tag.split('}')[-1]
g={};p={};calls=[];results=[]
def field(b,n):return b.find("field[@name='%s']"%n).text or ''
def value(b,n,env):return run(b.find("value[@name='%s']/block"%n),env)
def body(b,n,env):return run(b.find("statement[@name='%s']/block"%n),env)
def run(b,env):
 if b is None:return None
 t=b.attrib['type'];m=b.find('mutation');a=m.attrib if m is not None else {}
 if t=='text':return field(b,'TEXT')
 if t=='math_number':return int(field(b,'NUM'))
 if t=='logic_boolean':return field(b,'BOOL')=='TRUE'
 if t=='lexical_variable_get':
  name=field(b,'VAR');return g[name[7:]] if name.startswith('global ') else env[name]
 if t=='text_join':return ''.join(str(run(v.find('block'),env)) for v in b.findall('value'))
 if t=='text_length':return len(value(b,'VALUE',env))
 if t=='text_starts_at':return value(b,'TEXT',env).find(value(b,'PIECE',env))+1
 if t=='text_segment':
  text=value(b,'TEXT',env);start=value(b,'START',env)-1;return text[start:start+value(b,'LENGTH',env)]
 if t=='math_subtract':return value(b,'A',env)-value(b,'B',env)
 if t=='math_add':return sum(run(v.find('block'),env) for v in b.findall('value'))
 if t in ('logic_compare','math_compare'):
  x,y=value(b,'A',env),value(b,'B',env);op=field(b,'OP');return {'EQ':lambda:x==y,'GT':lambda:x>y}[op]()
 if t=='component_set_get':
  key=(a['instance_name'],a['property_name'])
  if a['set_or_get']=='get':return p.get(key)
  p[key]=value(b,'VALUE',env)
 elif t=='global_declaration':g[field(b,'NAME')]=value(b,'VALUE',env)
 elif t=='lexical_variable_set':g[field(b,'VAR')[7:]]=value(b,'VALUE',env)
 elif t=='component_method':calls.append((a['instance_name'],a['method_name'],[run(v.find('block'),env) for v in b.findall('value')]))
 elif t=='controls_if':body(b,'DO0' if value(b,'IF0',env) else 'ELSE',env)
 elif t=='controls_forEach':
  for item in value(b,'LIST',env):body(b,'DO',dict(env,**{field(b,'VAR'):item}))
 elif t=='component_event':body(b,'DO',env)
 else:raise AssertionError('Unknown '+t)
 run(b.find('next/block'),env)
events={}
for b in r.findall('block'):
 if b.attrib['type']=='global_declaration':run(b,{})
 if b.attrib['type']=='component_event':
  m=b.find('mutation').attrib;key=(m['instance_name'],m['event_name']);assert key not in events;events[key]=b
ble=json.loads(z.read('assets/external_comps/edu.mit.appinventor.ble/components.json'))[0]
methods={x['name']:x for x in ble['methods']};evs={x['name']:x for x in ble['events']}
for b in r.iter('block'):
 m=b.find('mutation')
 if m is not None and m.attrib.get('component_type')=='BluetoothLE':
  if b.attrib['type']=='component_method':assert len(b.findall('value'))==len(methods[m.attrib['method_name']]['params'])
  if b.attrib['type']=='component_event':assert m.attrib['event_name'] in evs
results.append('BLE method names, event names and argument counts match bundled extension metadata')
def fire(instance,name,**env):run(events[(instance,name)],env)
fire('Screen1','Initialize');assert p['BtnDisconnect','Enabled']==False
fire('BtnScan','Click');assert p['ScanClock','TimerEnabled']==True
fire('BluetoothLE1','DeviceFound');assert ('ListDevices','ElementsFromString') in p
p['ListDevices','SelectionIndex']=7;fire('ListDevices','AfterPicking');assert calls[-1]==('BluetoothLE1','Connect',[7]);assert p['ListDevices','Visible']==False
fire('BluetoothLE1','Connected');assert g['Session']==1;assert calls[-1][1]=='RegisterForStrings';assert calls[-1][2][-1]==False
results.append('initialize, scan, device selection and notification subscription')
p['BluetoothLE1','IsDeviceConnected']=True
fire('BluetoothLE1','StringsReceived',stringValues=['*Q12*','*Q13*']);assert g['Rx']=='*Q12**Q13*'
fire('WebViewer1','WebViewStringChange',value='POLL');assert p['WebViewer1','WebViewString']=='1|0|已連線|*Q12**Q13*';assert g['Rx']==''
fire('WebViewer1','WebViewStringChange',value='POLL');assert p['WebViewer1','WebViewString']=='1|0|已連線|'
results.append('multiple chunks accumulate and drain once on explicit webview poll')
fire('BluetoothLE1','StringsReceived',stringValues=['x'*32769]);assert g['Overflow']==1;assert g['Rx']==''
results.append('native bounded buffer overflow is counted and signaled')
csv='received_sample_index,adc_count\n1,12\n';fire('WebViewer1','WebViewStringChange',value='CSV:'+csv);assert calls[-1]==('File1','SaveFile',[csv,'Group2_ADC_latest100.csv']);assert p['WebViewer1','WebViewString']=='ACK'
fire('File1','AfterFileSaved',fileName='Group2_ADC_latest100.csv');assert calls[-1]==('Sharing1','ShareFile',['Group2_ADC_latest100.csv'])
results.append('CSV command strips prefix, writes correct text then opens sharing chooser')
fire('BluetoothLE1','Disconnected');assert p['BtnScan','Enabled']==True;assert g['Status']=='已中斷連線';assert p['ListDevices','Visible']==True
fire('BluetoothLE1','ConnectionFailed',reason='timeout');assert 'timeout' in calls[-1][2][0]
results.append('disconnect and connection failure restore controls')
(R/'native_block_test_results.json').write_text(json.dumps({'checks':results,'status':'passed','environment':'XML control-flow interpreter with mock components; not MIT import/build or Android hardware'},ensure_ascii=False,indent=2));print(len(results),'native control-flow checks passed')
