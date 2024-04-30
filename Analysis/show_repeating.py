#!/usr/bin/python3

import os
import sys
import time
import argparse
import subprocess

COMMAND = ""

def runCommand(cmd, args, file):
    process = subprocess.Popen(
        [cmd] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    command_in_loop = 0
    success = True
    while True:
        # if command_in_loop > 1e6:
        #     return 1
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            if "Found repeating pattern" in output or "[Error]" in output:
                print(output.strip())
                file.write(output.strip())
                file.write("\n")
            # command_in_loop += 1
            if "[Error]" in output:
                success = False

    exit_code = process.wait() if success else 1
    return exit_code

def runRepeating(cmd, args, file):
    count = 0
    while True:
        print()
        exit_code = runCommand(cmd, args, file)
        # time.sleep(1)
        if exit_code != 0:
            print(f"\nExit code {exit_code} returned")
            break
        count += 1
        os.system('clear')
        print(f"\rRun: {count} successful", end='', flush=True)
        file.write(f"Run: {count} successful\n")
        file.flush()
        # time.sleep(1)

    print(f"Completed runs: {count}")
    file.write(f"Completed runs: {count}\n")

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

    out_file = open("./Analysis/repeating.txt", "w")
    print()
    print()
    runRepeating(command, command_args, out_file)
    out_file.close()

if __name__ == "__main__":
    main()
