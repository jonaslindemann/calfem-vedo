# -*- coding: iso-8859-15 -*-
'''
ObjectiveFrame generated CALFEM Python code.
'''

import numpy as np
import calfem.core as cfc


# ---- Embedded df3 file

df3_file = '''MQoxCmRlZmF1bHQKNTIKMCAwIC0xCjIuMWUrMDkgNy45M2UrMTAgMC4xMjU2NjQgMC4wMDEyNTY2
NCAwLjAwMTI1NjY0IDAuMDAyNTEzMjcgCjYKMTEKMC4xIDAuMSAwLjEgMC4xIDAuMDEgMC4wMSAw
LjAxIDAuMSAwLjEgMC4yIDAuMDkgCjYKMi4xZSswOSAwLjEyNTY2NCA3LjkzZSsxMCAwLjAwMTI1
NjY0IDAuMDAxMjU2NjQgMC4wMDI1MTMyNyAKMTMKMCAwLjIgMAowLjEgMC4xNzMyMDUgMAowLjE3
MzIwNSAwLjEgMAowLjIgMCAwCjAuMTczMjA1IC0wLjEgMAowLjEgLTAuMTczMjA1IDAKMCAtMC4y
IDAKLTAuMSAtMC4xNzMyMDUgMAotMC4xNzMyMDUgLTAuMSAwCi0wLjIgMCAwCi0wLjE3MzIwNSAw
LjEgMAotMC4xIDAuMTczMjA1IDAKMCAwLjIgMAoKCjMyCjEgMSAtMC42IDAgLTIuNiAKMSAyIDAu
NCAwIC0yLjYgCjEgMyAwLjQgMCAtMS42IAoxIDQgLTAuNiAwIC0xLjYgCjEgNSAtMC42IDEgLTIu
NiAKMSA2IC0wLjYgMSAtMS42IAoxIDcgMC40IDEgLTIuNiAKMSA4IC0wLjYgMiAtMi42IAoxIDkg
LTAuNiAyIC0xLjYgCjEgMTAgMC40IDIgLTIuNiAKMSAxMSAwLjQgMSAtMS42IAoxIDEyIDAuNCAy
IC0xLjYgCjEgMTMgLTAuNiAxIC0wLjYgCjEgMTQgLTAuNiAyIC0wLjYgCjEgMTUgMC40IDIgLTAu
NiAKMSAxNiAwLjQgMSAtMC42IAoxIDE3IC0wLjYgMSAwLjQgCjEgMTggLTAuNiAyIDAuNCAKMSAx
OSAwLjQgMiAwLjQgCjEgMjAgMC40IDEgMC40IAoxIDIxIC0wLjYgMiAxLjQgCjEgMjIgMC40IDIg
MS40IAoxIDIzIC0wLjYgMSAxLjQgCjEgMjQgMC40IDEgMS40IAoxIDI1IC0wLjYgMiAyLjQgCjEg
MjYgMC40IDIgMi40IAoxIDI3IC0wLjYgMSAyLjQgCjEgMjggMC40IDEgMi40IAoxIDI5IC0wLjYg
MCAxLjQgCjEgMzAgMC40IDAgMS40IAoxIDMxIDAuNCAwIDIuNCAKMSAzMiAtMC42IDAgMi40IAoK
Cjg5CjEgIDIgMSA1IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEK
MiAgMiA1IDggCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQozICAy
IDIgNyAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjQgIDIgNyAx
MCAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjUgIDIgNCA2IAow
IDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNiAgMiA2IDkgCjAgMCAw
CjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo3ICAyIDMgMTEgCjAgMCAwCjAK
CjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo4ICAyIDExIDEyIAowIDAgMAowCgox
CjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKOSAgMiAxMyAxNCAKMCAwIDAKMAoKMQow
IDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjEwICAyIDE2IDE1IAowIDAgMAowCgoxCjAg
MCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMTEgIDIgMTcgMTggCjAgMCAwCjAKCjEKMCAw
IDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQoxMiAgMiAyMCAxOSAKMCAwIDAKMAoKMQowIDAg
MCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjEzICAyIDIzIDIxIAowIDAgMAowCgoxCjAgMCAw
IDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMTQgIDIgMjQgMjIgCjAgMCAwCjAKCjEKMCAwIDAg
MCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQoxNSAgMiAyOSAyMyAKMCAwIDAKMAoKMQowIDAgMCAw
IDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjE2ICAyIDMwIDI0IAowIDAgMAowCgoxCjAgMCAwIDAg
MCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMTcgIDIgMzIgMjcgCjAgMCAwCjAKCjEKMCAwIDAgMCAw
IDAgCjEKMCAwIDAgMCAwIDAgCjAgMQoxOCAgMiAzMSAyOCAKMCAwIDAKMAoKMQowIDAgMCAwIDAg
MCAKMQowIDAgMCAwIDAgMCAKMCAxCjE5ICAyIDI3IDI1IAowIDAgMAowCgoxCjAgMCAwIDAgMCAw
IAoxCjAgMCAwIDAgMCAwIAowIDEKMjAgIDIgMjggMjYgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAg
CjEKMCAwIDAgMCAwIDAgCjAgMQoyMSAgMiA1IDYgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEK
MCAwIDAgMCAwIDAgCjAgMQoyMiAgMiA3IDExIAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAg
MCAwIDAgMCAwIAowIDEKMjMgIDIgNiAxMSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAg
MCAwIDAgMCAKMCAxCjI0ICAyIDUgNyAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAw
IDAgMCAKMCAxCjI1ICAyIDggMTAgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAw
IDAgCjAgMQoyNiAgMiA5IDEyIAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAw
IAowIDEKMjcgIDIgOCA5IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAow
IDEKMjggIDIgMTAgMTIgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAg
MQoyOSAgMiA5IDE0IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEK
MzAgIDIgMTIgMTUgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQoz
MSAgMiA2IDEzIAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMzIg
IDIgMTYgMTMgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQozMyAg
MiAxMSAxNiAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjM0ICAy
IDE0IDE1IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMzUgIDIg
MTQgMTggCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQozNiAgMiAx
NSAxOSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjM3ICAyIDEz
IDE3IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKMzggIDIgMTYg
MjAgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQozOSAgMiAxOCAx
OSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjQwICAyIDE3IDIw
IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNDEgIDIgMTcgMjMg
CjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo0MiAgMiAyMCAyNCAK
MCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjQzICAyIDE5IDIyIAow
IDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNDQgIDIgMTggMjEgCjAg
MCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo0NSAgMiAyMSAyNSAKMCAw
IDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjQ2ICAyIDIyIDI2IAowIDAg
MAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNDcgIDIgMjQgMjggCjAgMCAw
CjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo0OCAgMiAyMyAyNyAKMCAwIDAK
MAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjQ5ICAyIDI4IDI3IAowIDAgMAow
CgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNTAgIDIgMjUgMjYgCjAgMCAwCjAK
CjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo1MSAgMiAyMiAyMSAKMCAwIDAKMAoK
MQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjUyICAyIDI0IDIzIAowIDAgMAowCgox
CjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNTMgIDIgMjYgMjQgCjAgMCAwCjAKCjEK
MCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo1NCAgMiAyMiAyMCAKMCAwIDAKMAoKMQow
IDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjU1ICAyIDE5IDE2IAowIDAgMAowCgoxCjAg
MCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNTYgIDIgMjAgMTUgCjAgMCAwCjAKCjEKMCAw
IDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo1NyAgMiAxNiAxMiAKMCAwIDAKMAoKMQowIDAg
MCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjU4ICAyIDExIDEwIAowIDAgMAowCgoxCjAgMCAw
IDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNTkgIDIgMzEgMjQgCjAgMCAwCjAKCjEKMCAwIDAg
MCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo2MCAgMiAyOCAzMCAKMCAwIDAKMAoKMQowIDAgMCAw
IDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjYxICAyIDMyIDI4IAowIDAgMAowCgoxCjAgMCAwIDAg
MCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKNjIgIDIgMjcgMzEgCjAgMCAwCjAKCjEKMCAwIDAgMCAw
IDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo2MyAgMiAyMyAzMiAKMCAwIDAKMAoKMQowIDAgMCAwIDAg
MCAKMQowIDAgMCAwIDAgMCAKMCAxCjY0ICAyIDI5IDI3IAowIDAgMAowCgoxCjAgMCAwIDAgMCAw
IAoxCjAgMCAwIDAgMCAwIAowIDEKNjUgIDIgMjcgMjYgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAg
CjEKMCAwIDAgMCAwIDAgCjAgMQo2NiAgMiAyOCAyNSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAK
MQowIDAgMCAwIDAgMCAKMCAxCjY3ICAyIDIzIDI1IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAox
CjAgMCAwIDAgMCAwIAowIDEKNjggIDIgMjYgMjEgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEK
MCAwIDAgMCAwIDAgCjAgMQo2OSAgMiAxNyAyMSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQow
IDAgMCAwIDAgMCAKMCAxCjcwICAyIDEzIDE4IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAg
MCAwIDAgMCAwIAowIDEKNzEgIDIgMTQgMTcgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAw
IDAgMCAwIDAgCjAgMQo3MiAgMiA5IDEzIAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAw
IDAgMCAwIAowIDEKNzMgIDIgOCA2IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAg
MCAwIAowIDEKNzQgIDIgMTAgNSAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAg
MCAKMCAxCjc1ICAyIDggMTIgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAg
CjAgMQo3NiAgMiA5IDE1IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAow
IDEKNzcgIDIgMTQgMTkgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAg
MQo3OCAgMiAxOCAyMiAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAx
Cjc5ICAyIDUgMTEgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo4
MCAgMiA2IDE2IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKODEg
IDIgMTMgMjAgCjAgMCAwCjAKCjEKMCAwIDAgMCAwIDAgCjEKMCAwIDAgMCAwIDAgCjAgMQo4MiAg
MiAxNyAyNCAKMCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjgzICAy
IDIzIDI4IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKODQgIDIg
MSA2IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKODUgIDIgMyA3
IAowIDAgMAowCgoxCjAgMCAwIDAgMCAwIAoxCjAgMCAwIDAgMCAwIAowIDEKODYgIDIgMTEgMiAK
MCAwIDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjg3ICAyIDcgMSAKMCAw
IDAKMAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjg4ICAyIDIgNSAKMCAwIDAK
MAoKMQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCjg5ICAyIDUgNCAKMCAwIDAKMAoK
MQowIDAgMCAwIDAgMCAKMQowIDAgMCAwIDAgMCAKMCAxCgoKMQoxIDEKNAoxOQoxNQoxNAoxOAow
IC0xMDAwIDAgCnNpbmdsZQoxCgoKMAoKCjIKOAoxCjIKMwo0CjI5CjMwCjMxCjMyCgoxIDAKMSAw
CjEgMAoxIDAKMSAwCjEgMApmaXhlZCBwb3Mvcm90CjEKMAoKMSAwCjEgMAoxIDAKMCAwCjAgMAow
IDAKZml4ZWQgcG9zCjEK'''

# ---- Materials

ep = np.array([
	[2.1e+09, 7.93e+10, 0.125664, 0.00125664, 0.00125664, 0.00251327],
])

# ---- Beam material idx

beam_mat_idx = np.array(
	[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
)


# ---- Bar material idx

bar_mat_idx = np.array(
	[],
)


# ---- Beam Edof

edof_beams = np.array([
	[1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 30],
	[25, 26, 27, 28, 29, 30, 43, 44, 45, 46, 47, 48],
	[7, 8, 9, 10, 11, 12, 37, 38, 39, 40, 41, 42],
	[37, 38, 39, 40, 41, 42, 55, 56, 57, 58, 59, 60],
	[19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36],
	[31, 32, 33, 34, 35, 36, 49, 50, 51, 52, 53, 54],
	[13, 14, 15, 16, 17, 18, 61, 62, 63, 64, 65, 66],
	[61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72],
	[73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84],
	[91, 92, 93, 94, 95, 96, 85, 86, 87, 88, 89, 90],
	[97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108],
	[115, 116, 117, 118, 119, 120, 109, 110, 111, 112, 113, 114],
	[133, 134, 135, 136, 137, 138, 121, 122, 123, 124, 125, 126],
	[139, 140, 141, 142, 143, 144, 127, 128, 129, 130, 131, 132],
	[169, 170, 171, 172, 173, 174, 133, 134, 135, 136, 137, 138],
	[175, 176, 177, 178, 179, 180, 139, 140, 141, 142, 143, 144],
	[187, 188, 189, 190, 191, 192, 157, 158, 159, 160, 161, 162],
	[181, 182, 183, 184, 185, 186, 163, 164, 165, 166, 167, 168],
	[157, 158, 159, 160, 161, 162, 145, 146, 147, 148, 149, 150],
	[163, 164, 165, 166, 167, 168, 151, 152, 153, 154, 155, 156],
	[25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
	[37, 38, 39, 40, 41, 42, 61, 62, 63, 64, 65, 66],
	[31, 32, 33, 34, 35, 36, 61, 62, 63, 64, 65, 66],
	[25, 26, 27, 28, 29, 30, 37, 38, 39, 40, 41, 42],
	[43, 44, 45, 46, 47, 48, 55, 56, 57, 58, 59, 60],
	[49, 50, 51, 52, 53, 54, 67, 68, 69, 70, 71, 72],
	[43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54],
	[55, 56, 57, 58, 59, 60, 67, 68, 69, 70, 71, 72],
	[49, 50, 51, 52, 53, 54, 79, 80, 81, 82, 83, 84],
	[67, 68, 69, 70, 71, 72, 85, 86, 87, 88, 89, 90],
	[31, 32, 33, 34, 35, 36, 73, 74, 75, 76, 77, 78],
	[91, 92, 93, 94, 95, 96, 73, 74, 75, 76, 77, 78],
	[61, 62, 63, 64, 65, 66, 91, 92, 93, 94, 95, 96],
	[79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90],
	[79, 80, 81, 82, 83, 84, 103, 104, 105, 106, 107, 108],
	[85, 86, 87, 88, 89, 90, 109, 110, 111, 112, 113, 114],
	[73, 74, 75, 76, 77, 78, 97, 98, 99, 100, 101, 102],
	[91, 92, 93, 94, 95, 96, 115, 116, 117, 118, 119, 120],
	[103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114],
	[97, 98, 99, 100, 101, 102, 115, 116, 117, 118, 119, 120],
	[97, 98, 99, 100, 101, 102, 133, 134, 135, 136, 137, 138],
	[115, 116, 117, 118, 119, 120, 139, 140, 141, 142, 143, 144],
	[109, 110, 111, 112, 113, 114, 127, 128, 129, 130, 131, 132],
	[103, 104, 105, 106, 107, 108, 121, 122, 123, 124, 125, 126],
	[121, 122, 123, 124, 125, 126, 145, 146, 147, 148, 149, 150],
	[127, 128, 129, 130, 131, 132, 151, 152, 153, 154, 155, 156],
	[139, 140, 141, 142, 143, 144, 163, 164, 165, 166, 167, 168],
	[133, 134, 135, 136, 137, 138, 157, 158, 159, 160, 161, 162],
	[163, 164, 165, 166, 167, 168, 157, 158, 159, 160, 161, 162],
	[145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156],
	[127, 128, 129, 130, 131, 132, 121, 122, 123, 124, 125, 126],
	[139, 140, 141, 142, 143, 144, 133, 134, 135, 136, 137, 138],
	[151, 152, 153, 154, 155, 156, 139, 140, 141, 142, 143, 144],
	[127, 128, 129, 130, 131, 132, 115, 116, 117, 118, 119, 120],
	[109, 110, 111, 112, 113, 114, 91, 92, 93, 94, 95, 96],
	[115, 116, 117, 118, 119, 120, 85, 86, 87, 88, 89, 90],
	[91, 92, 93, 94, 95, 96, 67, 68, 69, 70, 71, 72],
	[61, 62, 63, 64, 65, 66, 55, 56, 57, 58, 59, 60],
	[181, 182, 183, 184, 185, 186, 139, 140, 141, 142, 143, 144],
	[163, 164, 165, 166, 167, 168, 175, 176, 177, 178, 179, 180],
	[187, 188, 189, 190, 191, 192, 163, 164, 165, 166, 167, 168],
	[157, 158, 159, 160, 161, 162, 181, 182, 183, 184, 185, 186],
	[133, 134, 135, 136, 137, 138, 187, 188, 189, 190, 191, 192],
	[169, 170, 171, 172, 173, 174, 157, 158, 159, 160, 161, 162],
	[157, 158, 159, 160, 161, 162, 151, 152, 153, 154, 155, 156],
	[163, 164, 165, 166, 167, 168, 145, 146, 147, 148, 149, 150],
	[133, 134, 135, 136, 137, 138, 145, 146, 147, 148, 149, 150],
	[151, 152, 153, 154, 155, 156, 121, 122, 123, 124, 125, 126],
	[97, 98, 99, 100, 101, 102, 121, 122, 123, 124, 125, 126],
	[73, 74, 75, 76, 77, 78, 103, 104, 105, 106, 107, 108],
	[79, 80, 81, 82, 83, 84, 97, 98, 99, 100, 101, 102],
	[49, 50, 51, 52, 53, 54, 73, 74, 75, 76, 77, 78],
	[43, 44, 45, 46, 47, 48, 31, 32, 33, 34, 35, 36],
	[55, 56, 57, 58, 59, 60, 25, 26, 27, 28, 29, 30],
	[43, 44, 45, 46, 47, 48, 67, 68, 69, 70, 71, 72],
	[49, 50, 51, 52, 53, 54, 85, 86, 87, 88, 89, 90],
	[79, 80, 81, 82, 83, 84, 109, 110, 111, 112, 113, 114],
	[103, 104, 105, 106, 107, 108, 127, 128, 129, 130, 131, 132],
	[25, 26, 27, 28, 29, 30, 61, 62, 63, 64, 65, 66],
	[31, 32, 33, 34, 35, 36, 91, 92, 93, 94, 95, 96],
	[73, 74, 75, 76, 77, 78, 115, 116, 117, 118, 119, 120],
	[97, 98, 99, 100, 101, 102, 139, 140, 141, 142, 143, 144],
	[133, 134, 135, 136, 137, 138, 163, 164, 165, 166, 167, 168],
	[1, 2, 3, 4, 5, 6, 31, 32, 33, 34, 35, 36],
	[13, 14, 15, 16, 17, 18, 37, 38, 39, 40, 41, 42],
	[61, 62, 63, 64, 65, 66, 7, 8, 9, 10, 11, 12],
	[37, 38, 39, 40, 41, 42, 1, 2, 3, 4, 5, 6],
	[7, 8, 9, 10, 11, 12, 25, 26, 27, 28, 29, 30],
	[25, 26, 27, 28, 29, 30, 19, 20, 21, 22, 23, 24],
])


# ---- Bar Edof

edof_bars = np.array([
])


# ----- Nodes

coords_beams = np.array([
	[-0.6, -2.6, 0],
	[0.4, -2.6, 0],
	[0.4, -1.6, 0],
	[-0.6, -1.6, 0],
	[-0.6, -2.6, 1],
	[-0.6, -1.6, 1],
	[0.4, -2.6, 1],
	[-0.6, -2.6, 2],
	[-0.6, -1.6, 2],
	[0.4, -2.6, 2],
	[0.4, -1.6, 1],
	[0.4, -1.6, 2],
	[-0.6, -0.6, 1],
	[-0.6, -0.6, 2],
	[0.4, -0.6, 2],
	[0.4, -0.6, 1],
	[-0.6, 0.4, 1],
	[-0.6, 0.4, 2],
	[0.4, 0.4, 2],
	[0.4, 0.4, 1],
	[-0.6, 1.4, 2],
	[0.4, 1.4, 2],
	[-0.6, 1.4, 1],
	[0.4, 1.4, 1],
	[-0.6, 2.4, 2],
	[0.4, 2.4, 2],
	[-0.6, 2.4, 1],
	[0.4, 2.4, 1],
	[-0.6, 1.4, 0],
	[0.4, 1.4, 0],
	[0.4, 2.4, 0],
	[-0.6, 2.4, 0],
])

coords_bars = np.array([
	[-0.6, -2.6, 0],
	[0.4, -2.6, 0],
	[0.4, -1.6, 0],
	[-0.6, -1.6, 0],
	[-0.6, -2.6, 1],
	[-0.6, -1.6, 1],
	[0.4, -2.6, 1],
	[-0.6, -2.6, 2],
	[-0.6, -1.6, 2],
	[0.4, -2.6, 2],
	[0.4, -1.6, 1],
	[0.4, -1.6, 2],
	[-0.6, -0.6, 1],
	[-0.6, -0.6, 2],
	[0.4, -0.6, 2],
	[0.4, -0.6, 1],
	[-0.6, 0.4, 1],
	[-0.6, 0.4, 2],
	[0.4, 0.4, 2],
	[0.4, 0.4, 1],
	[-0.6, 1.4, 2],
	[0.4, 1.4, 2],
	[-0.6, 1.4, 1],
	[0.4, 1.4, 1],
	[-0.6, 2.4, 2],
	[0.4, 2.4, 2],
	[-0.6, 2.4, 1],
	[0.4, 2.4, 1],
	[-0.6, 1.4, 0],
	[0.4, 1.4, 0],
	[0.4, 2.4, 0],
	[-0.6, 2.4, 0],
])

dofs_beams = np.array([
	[1, 2, 3, 4, 5, 6],
	[7, 8, 9, 10, 11, 12],
	[13, 14, 15, 16, 17, 18],
	[19, 20, 21, 22, 23, 24],
	[25, 26, 27, 28, 29, 30],
	[31, 32, 33, 34, 35, 36],
	[37, 38, 39, 40, 41, 42],
	[43, 44, 45, 46, 47, 48],
	[49, 50, 51, 52, 53, 54],
	[55, 56, 57, 58, 59, 60],
	[61, 62, 63, 64, 65, 66],
	[67, 68, 69, 70, 71, 72],
	[73, 74, 75, 76, 77, 78],
	[79, 80, 81, 82, 83, 84],
	[85, 86, 87, 88, 89, 90],
	[91, 92, 93, 94, 95, 96],
	[97, 98, 99, 100, 101, 102],
	[103, 104, 105, 106, 107, 108],
	[109, 110, 111, 112, 113, 114],
	[115, 116, 117, 118, 119, 120],
	[121, 122, 123, 124, 125, 126],
	[127, 128, 129, 130, 131, 132],
	[133, 134, 135, 136, 137, 138],
	[139, 140, 141, 142, 143, 144],
	[145, 146, 147, 148, 149, 150],
	[151, 152, 153, 154, 155, 156],
	[157, 158, 159, 160, 161, 162],
	[163, 164, 165, 166, 167, 168],
	[169, 170, 171, 172, 173, 174],
	[175, 176, 177, 178, 179, 180],
	[181, 182, 183, 184, 185, 186],
	[187, 188, 189, 190, 191, 192],
])

dofs_bars = np.array([
	[1, 2, 3],
	[7, 8, 9],
	[13, 14, 15],
	[19, 20, 21],
	[25, 26, 27],
	[31, 32, 33],
	[37, 38, 39],
	[43, 44, 45],
	[49, 50, 51],
	[55, 56, 57],
	[61, 62, 63],
	[67, 68, 69],
	[73, 74, 75],
	[79, 80, 81],
	[85, 86, 87],
	[91, 92, 93],
	[97, 98, 99],
	[103, 104, 105],
	[109, 110, 111],
	[115, 116, 117],
	[121, 122, 123],
	[127, 128, 129],
	[133, 134, 135],
	[139, 140, 141],
	[145, 146, 147],
	[151, 152, 153],
	[157, 158, 159],
	[163, 164, 165],
	[169, 170, 171],
	[175, 176, 177],
	[181, 182, 183],
	[187, 188, 189],
])


# ----- Local Z-axis

beam_orientation = np.array([
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[-6.12323e-17, 1, 1.22465e-16],
	[1, 0, 0],
	[1, 0, 0],
	[6.12323e-17, -1, -6.12323e-17],
	[6.12323e-17, -1, -6.12323e-17],
	[6.12323e-17, -1, -6.12323e-17],
	[6.12323e-17, -1, -6.12323e-17],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[6.12323e-17, -1, 6.12323e-17],
	[1, 0, 0],
	[6.12323e-17, -1, -6.12323e-17],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[6.12323e-17, -1, -6.12323e-17],
	[6.12323e-17, -1, -6.12323e-17],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[1, 0, 0],
	[6.12323e-17, -1, 6.12323e-17],
	[6.12323e-17, -1, -6.12323e-17],
	[6.12323e-17, -1, 6.12323e-17],
	[6.12323e-17, -1, 6.12323e-17],
	[-1.47828e-16, -0.707107, 0.707107],
	[-1.47828e-16, -0.707107, 0.707107],
	[-1.47828e-16, -0.707107, 0.707107],
	[6.12323e-17, -0.707107, -0.707107],
	[6.12323e-17, -0.707107, -0.707107],
	[6.12323e-17, -0.707107, -0.707107],
	[6.12323e-17, -0.707107, -0.707107],
	[-1.47828e-16, -0.707107, 0.707107],
	[5.55112e-17, -1, -5.55112e-17],
	[0, -1, -1.11022e-16],
	[2.53633e-17, -0.707107, -0.707107],
	[6.12323e-17, -0.707107, 0.707107],
	[5.55112e-17, -1, -5.55112e-17],
	[5.55112e-17, -1, 5.55112e-17],
	[6.12323e-17, -0.707107, 0.707107],
	[0.707107, -0.707107, 1.79345e-17],
	[6.12323e-17, -0.707107, 0.707107],
	[6.12323e-17, -0.707107, 0.707107],
	[2.53633e-17, -0.707107, -0.707107],
	[2.53633e-17, -0.707107, -0.707107],
	[2.53633e-17, -0.707107, -0.707107],
	[-1.66533e-16, -1, -5.55112e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[0.707107, -0.707107, -1.79345e-17],
	[6.12323e-17, -0.707107, 0.707107],
	[6.12323e-17, -0.707107, -0.707107],
	[-1.47828e-16, -0.707107, 0.707107],
	[-1.66533e-16, -1, -5.55112e-17],
	[5.55112e-17, -1, 5.55112e-17],
	[2.53633e-17, -0.707107, -0.707107],
])

# ----- Boundary conditions

node_bc = np.array([
	[1, 1, 1, 1, 1, 1],
	[1, 1, 1, 0, 0, 0],
])

node_bc_val = np.array([
	[0, 0, 0, 0, 0, 0],
	[0, 0, 0, 0, 0, 0],
])

node_bc_idx = np.array([
	[1, 0],
	[2, 0],
	[3, 0],
	[4, 0],
	[29, 0],
	[30, 0],
	[31, 0],
	[32, 0],
])

node_load = np.array([
	[0, 0, -1000],
])

node_load_idx = np.array([
	[14, 0],
	[15, 0],
	[18, 0],
	[19, 0],
])

