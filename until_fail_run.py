#!/usr/bin/python3

import os
import sys
import time
import argparse
import subprocess

COMMAND = ""

def runUntilFailure(cmd, args):
    count = 0
    while True:
        process = subprocess.Popen(
            [cmd] + args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        stdout, stderr = process.communicate()

        # exit_code = p.wait()
        exit_code = process.returncode

        if exit_code != 0:
            print(f"command: {cmd} failed.")
            print("Exit code: ", process.returncode)
            print("Standard Output:\n", stdout.decode())
            print("Error Output:\n", stderr.decode())
            break
        else:
            count += 1
            # print(f"Run {count} successful")
            print(f"\rRun: {count} successful", end='', flush=True)
            time.sleep(1)

    print(f"Completed runs: {count}")

def main():
    # parser = argparse.ArgumentParser(description="Run a command repeatedly until it fails.")

    # parser.add_argument("command_str",
    #                     type=str,
    #                     help="The command that should be run repeatedly.")
    # parser.add_argument("command_args", nargs='*', help="Additional arguments to pass to the command.")

    # args = parser.parse_args()

    # command = args.command_str
    # command_args = args.command_args

    # if not command:
    #     print(f"[Error]: Not command argument provided.\nExiting")
    #     exit(1)

    if len(sys.argv) < 2:
        print(f"[Error]: Not command argument provided.\nExiting")
        sys.exit(1)

    command = sys.argv[1]

    command_args = sys.argv[2:]

    if command[0:2] == "./":
        command = os.path.abspath(command)

    print()
    print()
    runUntilFailure(command, command_args)

if __name__ == "__main__":
    main()
