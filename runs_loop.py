#!/usr/bin/python3

import os
import sys
import time
import subprocess

COMMAND = ""

def runCommand(cmd, args):
    process = subprocess.Popen(
        [cmd] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )

    # stdout, stderr = process.communicate()

    # # exit_code = p.wait()
    # exit_code = process.returncode

    # # print("Standard Output:\n", stdout.decode())
    # # print("Error Output:\n", stderr.decode())
    # # print("Exit code: ", exit_code)
    # time.sleep(1)
    # os.system('clear')

    for line in process.stdout:
        if line:
            line = str(line.strip(), 'utf-8')
            if "[Error]" in line:
                sys.exit(1)
            print(line)
    time.sleep(1)
    os.system('clear')

    exit_code = 0
    return exit_code

def main():

    if len(sys.argv) < 2:
        print(f"[Error]: Not command argument provided.\nExiting")
        sys.exit(1)

    command = sys.argv[1]

    command_args = sys.argv[2:]

    if command[0:2] == "./":
        command = os.path.abspath(command)

    print()
    print()
    count = 0
    while True:
        if runCommand(command, command_args) != 0:
            break;
        time.sleep(1)
        count += 1
        print(f"Running #{count} iteration")

if __name__ == "__main__":
    main()

