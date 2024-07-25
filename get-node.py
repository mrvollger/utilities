#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import subprocess
import math


def run_cmd_and_get_stdout(cmd):
    result = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
    )
    return result.stdout.decode("utf-8").strip()


def get_allowed_time():
    """_summary_
    # get time till next maitnence
    NEXTRES=$(" )
    NOW=$(date +%s)
    TIME_UNTIL=$(((NEXTRES - NOW) / 3600))
    echo "$TIME_UNTIL hours left until reservation begins"
    TIME_REQ=$((TIME_UNTIL-1))
    echo "Requesting $TIME_REQ hours"
    """
    # find the time until the next reservation
    nextres = run_cmd_and_get_stdout(
        "scontrol show res | grep Maintenance | head -n1 | awk '{print $2}' | cut -f2 -d= | xargs -I {} date +%s --date \"{}\""
    )
    cur_time = run_cmd_and_get_stdout("date +%s")
    try:
        time_until = (int(nextres) - int(cur_time)) / 3600
        if time_until < 0:
            nextres = run_cmd_and_get_stdout(
                "scontrol show res | grep ReservationName | grep Maintenance | head -n2 | tail -n 1 | awk '{print $2}' | cut -f2 -d= | xargs -I {} date +%s --date \"{}\""
            )
            time_until = (int(nextres) - int(cur_time)) / 3600
        logging.info(f"next res minus cur time: {nextres/3600} - {cur_time/3600}")
    except:
        logging.error(
            f"Could not get time until next reservation: {nextres} - {cur_time}"
        )
    logging.info(f"{time_until:.2f} hours left until reservation begins")
    time_req = math.floor(time_until - 1 / 60)
    logging.info(f"Requesting {time_req} hours")
    return time_req


def get_job(node, partition, cores, mem, account, time, odir):
    if time is None:
        time = get_allowed_time()

    if node is None:
        node = ""
    else:
        node = f"-w {node}"

    run_cmd_and_get_stdout(f"mkdir -p {odir}")
    run_cmd_and_get_stdout(
        f"printf '#!/bin/bash\nsleep {time+1}h\n' > {odir}/tmp_node_script.sh"
    )
    comment = "'requested by get-node.py'"

    cmd = f"sbatch --output {odir}/tmp_node_script.out --parsable -A {account} --mem={mem} --time={time}:00:00 -p {partition} -c {cores} {odir}/tmp_node_script.sh  --comment={comment}"

    logging.info(f"sbatch request:\n{cmd}")
    job_id = run_cmd_and_get_stdout(cmd)
    logging.info(f"Job ID: {job_id}")


def main(
    *,
    node: Optional[str] = None,
    partition: str = "compute,cpu-g2,cpu-g2-mem2x",
    cores: int = 4,
    mem: str = "100G",
    account: str = "stergachislab",
    odir: Path = Path("~/get-node-output"),
    time: Optional[str] = None,
    verbose: int = 1,
):
    """
    Author Mitchell R. Vollger

    :param cores: The number of cores to request
    :param mem: The amount of memory to request
    :param time: The amount of time to request, if not provided will use the time until the next maintenance window.
    :param partition: The partition to request
    :param account: The account to use
    :param node: The node to request
    :param odir: The output directory to write the output to
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    get_job(node, partition, cores, mem, account, time, odir)

    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")


exit

# DIR=${USER}/get-node-output
# mkdir -p $DIR
# printf "#!/bin/bash\nsleep 360d\n" > ${DIR}/tmp_node_script.sh
# job_id=$()
