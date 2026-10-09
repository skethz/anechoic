#!/usr/bin/env python3
"""Record the built image and exact inputs; this never programs hardware."""
import hashlib, json, math, os, re, struct
from pathlib import Path

root=Path(os.environ['TASK_ROOT'])
board=root/'build'/os.environ.get('BOARD_TAG','board_'+os.environ['KERNEL_MHZ']+'mhz')
sfx=('_'+os.environ['SCA_VARIANT']) if os.environ.get('SCA_VARIANT') else ''
iprepo=Path(os.environ.get('SNOWBALL_IP_REPO',str(root/'build/iprepo')))
clock=json.loads((board/'kernel_clock.json').read_text())
assert math.isfinite(clock['actual_mhz']) and 0<clock['actual_mhz']<1000
assert abs(1000/clock['actual_mhz']-clock['routed_period_ns'])<.002
assert math.isclose(clock['actual_mhz'],clock['input_mhz']*clock['feedback_multiplier']/clock['input_divider']/clock['output_divider'],rel_tol=1e-8)
uuid_path=board/'prj.runs/impl_1/pfm_uuid_manifest.dict'
manifest=uuid_path.read_text()
uuid=re.search(r'\blogic_uuid\s+([0-9a-f]{32})\b',manifest).group(1)
checkpoint=board/'prj.runs/synth_1/top_wrapper.dcp'
assert hashlib.md5(checkpoint.read_bytes()).hexdigest()==uuid
elf=(board/'amc.elf').read_bytes()
assert elf[:6]==b'\x7fELF\x01\x01' and struct.unpack_from('<H',elf,18)[0]==40, 'Expected ARM ELF32 management firmware'
files={
 'pdi':board/'snowball_v80.pdi',
 'uuid_manifest':uuid_path,
 'hardware_pdi':board/'prj.runs/impl_1/top_wrapper.pdi',
 'management_firmware':board/'amc.elf',
 'hardware_handoff':board/'amd_v80_gen5x8_25.1.xsa',
 'routed_checkpoint':board/'prj.runs/impl_1/top_wrapper_routed.dcp',
 'routed_timing':board/'routed_timing.rpt',
 'timing_signoff':board/'timing_signoff.json',
 'timer_structure':board/'timer_structure.json',
 'bus_skew':board/'bus_skew.rpt',
 'drc':board/'drc.rpt',
 'host_binary':root/('build/host_sca'+sfx),
 'graph':root/'data/K2000.bin',
 'instrumented_ip':iprepo/'sca_v80/component.xml',
 'counter_rtl':iprepo/'sca_v80/hdl/verilog/sca_v80_control_s_axi.v',
 'instrumentation':root/('results/instrumentation'+sfx+'.json'),
}
records={}
for key,path in files.items():
 data=path.read_bytes()
 assert data, ('Empty build artifact',str(path))
 records[key]={'path':str(path),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
assert records['pdi']['bytes']<132120576, 'Image exceeds inventoried flash-partition capacity'
skew=(board/'bus_skew.rpt').read_text()
skew_rows=re.findall(r'^\s*(?:Slow|Fast)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$',skew,re.M)
assert skew_rows and min(float(row[2]) for row in skew_rows)>=0, 'Bus-skew verification failed'
result={'status':'Built image; physical execution not yet established','logic_uuid':uuid,
 'physical_board_execution':False,'clock':clock,'files':records,
 'pdi_format':'Normal boot image with matching AMC; no flash-partition-table prefix',
 'recovery_image':{'path':'/opt/amd/aved/amd_v80_gen5x8_25.1_exdes_1_xbtest_stress/design.pdi',
  'sha256':'03c299444401bb0ae98497134cbdccdbedd5a3d0f11709a6ede4c94eef248201',
  'logic_uuid':'5a73aca29b601abfd14ffca90557701b'}}
(board/'image_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
