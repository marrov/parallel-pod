# Query user for their 3DPOD variables
MSIZE = float(input("Input number of points: "))
TSIZE = float(input("Input number of snapshots: "))
VSIZE = float(input("Input number of components: "))
flag_spod = float(input("Select SPOD filtering (1) or normal POD (0): "))
flag_s_or_d = float(input("Select double (1) or single (0): "))

if flag_s_or_d == 1:
    s_or_d = 8
elif flag_s_or_d == 0:
    s_or_d = 4


m_size = (MSIZE*TSIZE*VSIZE*s_or_d)/(1024**3)
pm_size = (TSIZE*TSIZE*s_or_d)/(1024**3)
spm_size = ((TSIZE*TSIZE*s_or_d)/(1024**3))*flag_spod
pod_size = (MSIZE*TSIZE*VSIZE*s_or_d)/(1024**3)

print (f"""
RAM memory usage report for 3DPOD:

- m matrix occupies {round(m_size,2)} Gb of RAM
- pm matrix occupies {round(pm_size,2)} Gb of RAM
- spm matrix occupies {round(spm_size,2)} Gb of RAM
- pod matrix occupies {round(pod_size,2)} Gb of RAM

Total RAM usage is {round(m_size+pm_size+spm_size+pod_size,2)} Gb
""")
