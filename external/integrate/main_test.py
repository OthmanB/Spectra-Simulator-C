import numpy as np
from subprocess import Popen, PIPE

l=1
ftype='gauss'
theta0=np.pi/2
delta=np.pi/8

process = Popen(["./Alm", str(l), str(theta0), str(delta), ftype], stdout=PIPE, stderr=PIPE)
(output, err) = process.communicate()
exit_code = process.wait()


print("Output retrieved from python:")
print(output)
exit()

