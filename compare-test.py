def parse_liquid_tests(lines):
    """
    Example line:
         2[    0.04 ms]   PASS   passed    4/   4 checks (100.0%) agc_crcf_dc_gain_control

    Example module line:
    120: fskmodem:
    """
    import re
    from collections import defaultdict

    TEST_RE = re.compile(r"^\s*\d+\[\s*\d+\.\d+ ms\]\s+(\w+)\s+passed\s+\d+/\s*\d+\s+checks\s+\(\s*\d+\.\d+%\)\s+(\w+)")
    MODULE_RE = re.compile(r"^(\d+): (\w+):$")

    module_map = defaultdict(list)
    passing_tests = []
    failing_tests = []
    current_module = None

    for line in lines:
        if m := TEST_RE.match(line):
            test_name = m.group(2)
            module_map[current_module].append(test_name)
            if m.group(1) == "PASS":
                passing_tests.append(test_name)
            else:
                failing_tests.append(test_name)

        if m := MODULE_RE.match(line):
            current_module = m.group(2)

    print(len(module_map), len(passing_tests), len(failing_tests))

    return module_map, passing_tests, failing_tests


def parse_yagi_tests(lines):
    """
    Example annotation line:
        !! liquid test annotation: autotest_iirfiltsos_step_n2 -> test_iirfiltsos_step_n2 !!

    Example test line:
    running 2 tests
        test filter::iirfiltsos::tests::test_iirfiltsos_impulse_n2 ... ok
        test filter::iirfiltsos::tests::test_iirfiltsos_step_n2 ... FAILED
    """
    import re

    RUNNING_TEST_RE = re.compile(r"^running (\d+) tests$")
    TEST_RE = re.compile(r"^\s*test (\S+) \.\.\. (ok|FAILED)$")
    TEST_ANNOTATION_RE = re.compile(r"^\s*!! liquid test annotation: autotest_(\S+) -> (\S+) !!$")

    test_map = {}
    test_count = None
    passing_tests = []
    failing_tests = []

    for line in lines:
        if m := RUNNING_TEST_RE.match(line):
            test_count = int(m.group(1))

        if test_count:
            if m := TEST_RE.match(line):
                test_name = m.group(1).split("::")[-1]
                test_result = m.group(2)
                if test_result == "ok":
                    passing_tests.append(test_name)
                else:
                    failing_tests.append(test_name)
                test_count -= 1

        if m := TEST_ANNOTATION_RE.match(line):
            liquid_test_name = m.group(1)
            yagi_test_name = m.group(2)
            if yagi_test_name in test_map:
                print(f"Duplicate test name: {yagi_test_name}")
            test_map[yagi_test_name] = liquid_test_name

    liquid_passing = [test_map[test] for test in passing_tests if test in test_map]
    liquid_failing = [test_map[test] for test in failing_tests if test in test_map]

    return liquid_passing, liquid_failing


def main():
    import os

    with open("liquid-test.out", "r") as f:
        liquid_module_map, liquid_passing, liquid_failing = parse_liquid_tests(f.readlines())

    with open("yagi-test.out", "r") as f:
        yagi_passing, yagi_failing = parse_yagi_tests(f.readlines())

    all_liquid_pass = 0
    all_liquid_fail = 0
    all_liquid_missing = 0
    all_yagi_pass = 0
    all_yagi_fail = 0
    all_yagi_missing = 0

    for module, tests in liquid_module_map.items():
        liquid_pass = 0
        liquid_fail = 0
        liquid_missing = 0
        yagi_pass = 0
        yagi_fail = 0
        yagi_missing = 0
        for test in tests:
            if test in liquid_passing:
                liquid_pass += 1
            elif test in liquid_failing:
                liquid_fail += 1
            else:
                liquid_missing += 1

            if test in yagi_passing:
                yagi_pass += 1
            elif test in yagi_failing:
                yagi_fail += 1
            else:
                yagi_missing += 1

        all_liquid_pass += liquid_pass
        all_liquid_fail += liquid_fail
        all_liquid_missing += liquid_missing
        all_yagi_pass += yagi_pass
        all_yagi_fail += yagi_fail
        all_yagi_missing += yagi_missing

        yagi_pass_pct = yagi_pass / (liquid_pass + liquid_fail) * 100
        yagi_present_pct = (yagi_pass + yagi_fail) / (liquid_pass + liquid_fail) * 100

        print(f"{module:30s}: {liquid_pass:4d} {liquid_fail:4d} {liquid_missing:4d} {yagi_pass:4d} {yagi_fail:4d} {yagi_missing:4d} {yagi_pass_pct:8.2f}% {yagi_present_pct:8.2f}%")

    all_yagi_pass_pct = all_yagi_pass / (all_liquid_pass + all_liquid_fail) * 100
    all_yagi_present_pct = (all_yagi_pass + all_yagi_fail) / (all_liquid_pass + all_liquid_fail) * 100
    total = "TOTAL"
    print(f"\n{total:30s}: {all_liquid_pass:4d} {all_liquid_fail:4d} {all_liquid_missing:4d} {all_yagi_pass:4d} {all_yagi_fail:4d} {all_yagi_missing:4d} {all_yagi_pass_pct:8.2f}% {all_yagi_present_pct:8.2f}%")


    # write result to LIQUID_COMPAT.md
    with open("LIQUID_COMPAT.md", "w") as f:
        for module, tests in liquid_module_map.items():
            f.write(f"\n\n## {module}\n")
            f.write(f"| Test | Liquid | Yagi |\n")
            f.write(f"| ---- | ------ | ---- |\n")
            for test in tests:
                f.write(f"| {test:40s} |")
                if test in liquid_passing:
                    f.write(f" ✅")
                else:
                    f.write(f" ❌")

                f.write(f" |")
                if test in yagi_passing:
                    f.write(f" ✅")
                elif test in yagi_failing:
                    f.write(f" ❌")
                else:
                    f.write(f" ❓")
                f.write(f" |\n")

if __name__ == "__main__":
    main()
