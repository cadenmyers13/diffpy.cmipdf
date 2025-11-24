import argparse

from diffpy.cmipdf.version import __version__  # noqa


def main():
    parser = argparse.ArgumentParser(
        prog="diffpy.cmipdf",
        description=(
            "Python package for doing science.\n\n"
            "For more information, visit: "
            "https://github.com/diffpy/diffpy.cmipdf/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Show the program's version number and exit",
    )

    args = parser.parse_args()

    if args.version:
        print(f"diffpy.cmipdf {__version__}")
    else:
        # Default behavior when no arguments are given
        parser.print_help()


if __name__ == "__main__":
    main()
