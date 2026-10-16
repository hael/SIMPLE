"""Hardware discovery and usable calibration resources."""
import platform
from resource_planner.slurm import allocation
import os
from pathlib import Path

MIB = 1048576

def cpu_pool():
    sockets, seen = {}, set()
    for cpu in sorted(os.sched_getaffinity(0)):
        p = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        pair = (int((p/'physical_package_id').read_text()), int((p/'core_id').read_text()))
        if pair not in seen:
            sockets.setdefault(pair[0], []).append(cpu)
            seen.add(pair)
    if os.environ.get('SLURM_JOB_ID'):
        granted = int(os.environ.get('SLURM_CPUS_PER_TASK', '0'))
        if granted < 1:
            raise ValueError('Slurm calibration requires --cpus-per-task and a bound srun task')
        return [cpu for group in sockets.values() for cpu in group][:granted]
    # Reserve at least 20% of physical cores for desktop/OS use.
    return [cpu for group in sockets.values() for cpu in group[:max(1, int(len(group)*0.8))]]



def discover():
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    memory = {}
    if Path("/proc/meminfo").exists():
        for line in Path("/proc/meminfo").read_text().splitlines():
            key, value = line.split(":", 1)
            if key in ("MemTotal", "MemAvailable"):
                memory[key + "_bytes"] = int(value.split()[0]) * 1024
    return {"os": platform.system(), "architecture": platform.machine(),
            "logical_cpus": os.cpu_count(), "allowed_cpu_ids": affinity,
            "host_memory": memory, "allocation": allocation(),
            "warnings": ["Host inventory is not an allocation grant. Supply an explicit usable campaign budget.",
                         "Cgroup memory discovery, distributed scaling, NUMA and GPU qualification are not implemented."]}

