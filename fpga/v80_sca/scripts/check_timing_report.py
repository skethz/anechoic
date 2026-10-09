#!/usr/bin/env python3
"""Gate image packaging on the aggregate setup, hold and pulse-width report."""
import json,re,sys
from pathlib import Path
p=Path(sys.argv[1]);text=p.read_text();lines=text.splitlines()
header=next(i for i,s in enumerate(lines) if 'WNS(ns)' in s and 'WPWS(ns)' in s)
values=None
for line in lines[header+1:header+6]:
    fields=line.split()
    if len(fields)==12:
        try: values=[float(x) for x in fields];break
        except ValueError: pass
assert values is not None, 'Missing aggregate timing results'
assert all(values[i]>=0 for i in [0,4,8]), 'Setup, hold, or pulse-width timing failed'
assert all(values[i]==0 for i in [2,6,10]), 'Failing timing endpoints remain'
assert 'There are 0 pins that are not constrained for maximum delay.' in text
for check in ['no_clock','multiple_clock','generated_clocks','loops','latch_loops']:
    m=re.search(r'checking '+check+r' \((\d+)\)',text)
    assert m and int(m.group(1))==0, ('Invalid timing constraints',check)
result={'status':'PASS for constrained synchronous timing','setup_slack_ns':values[0],
 'hold_slack_ns':values[4],'pulse_width_slack_ns':values[8],
 'setup_failing_endpoints':int(values[2]),'hold_failing_endpoints':int(values[6]),
 'pulse_width_failing_endpoints':int(values[10]),'source_report':str(p),
 'external_delay_note':'External ports and constant-clock endpoints still require review of the verbose report; this is not a board-level protocol timing measurement.'}
p.with_name('timing_signoff.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
