"""Command Line Interface for using canuvit."""

import re
import click
from canuvit import observe_VIS, observe_UV, observe, __version__

INSTRUMENT_CHOICES = ("uvit", "sxt", "czti", "laxpc")
DEFAULT_INSTRUMENT = "uvit"

observe_func = {"all": observe, "vis": observe_VIS, "uv": observe_UV}

# Mutually exclusive click option: from https://stackoverflow.com/questions/49387833/prohibit-passing-several-feature-switches-in-python-click
class OnceSameNameOption(click.Option):
    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            if self.name in state.opts:
                param_same_name = [
                    opt.opts[0]
                    for opt in ctx.command.params
                    if isinstance(opt, OnceSameNameOption) and opt.name == self.name
                ]

                raise click.UsageError(
                    "Illegal usage: `{}` are mutually exclusive arguments.".format(
                        ", ".join(param_same_name)
                    )
                )

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OnceSameNameOption, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


class CheckREParamType(click.ParamType):
    def __init__(self, pattern: re.Pattern, name: str, humanpattern: str = None):
        super().__init__()
        self.name = name
        self.pattern = pattern
        self.humanpattern = humanpattern

    def convert(self, value, param, ctx):
        try:
            if self.pattern.match(value) is None:
                self.fail(
                    f"{value} is in an invalid format. Please use the {self.humanpattern} format.",
                    param,
                    ctx,
                )
        except TypeError:
            self.fail(f"{value} if of invalid type", param, ctx)
        return value


RAre = re.compile(r"^\d{2}:\d{2}:\d{2}(?:\.\d*)?$")
DECre = re.compile(r"^(?:\+|\-)?\d{2}:\d{2}:\d{2}(?:\.\d*)?$")

RAstr = CheckREParamType(RAre, "RA", "hh:mm:ss[.ss]")
DECstr = CheckREParamType(DECre, "DEC", "[-]dd:mm:ss[.ss]")


@click.command("canuvit", context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--all",
    "context",
    flag_value="all",
    default=True,
    help="Check safety for all filters. Default behaviour.",
    cls=OnceSameNameOption,
)
@click.option(
    "--vis",
    "context",
    flag_value="vis",
    help="Check saftey for only visible filters.",
    cls=OnceSameNameOption,
)
@click.option(
    "--uv",
    "context",
    flag_value="uv",
    help="Check safety for only UV filters.",
    cls=OnceSameNameOption,
)
@click.option(
    "-r",
    "--ra",
    metavar="RA",
    type=RAstr,
    required=True,
    help='Right ascension of the coordinate. Format: hh:mm:ss[.ss] e.g. "00:54:53.45"',
)
@click.option(
    "-d",
    "--dec",
    metavar="DEC",
    type=DECstr,
    required=True,
    help='Declination of the coordinate. Format: [-]dd:mm:ss[.ss] e.g. "-37:41:03.23".',
)
@click.option(
    "-i",
    "--instrument",
    type=click.Choice(INSTRUMENT_CHOICES),
    default=DEFAULT_INSTRUMENT,
    help=f"Instrument to check for.",
    show_default=True,
)
@click.option("-v", "--verbose", count=True, help="Increase output verbosity.")
@click.version_option(__version__, "--version", prog_name="canuvit")
def cli(context, ra, dec, instrument, verbose):
    """Program to check if a given coordinate can be safely observed using UVIT.

    \b
    Example usage:
    canuvit -r "13:12:14" -d "-14:15:13" """
    if verbose >= 1:
        click.echo(
            f"Checking {context} filters for coordinates: RA: {ra} DEC: {dec} in instrument: {instrument}."
        )
    observe_func[context](instrument, ra, dec)


if __name__ == "__main__":
    print("Please use either 'python -m canuvit' or just 'canuvit' to run the script.")
