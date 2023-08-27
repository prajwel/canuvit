#!/usr/bin/env python3


import warnings
import numpy as np
from io import BytesIO
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bs4 import BeautifulSoup
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
from matplotlib.colors import LogNorm
from requests import Session, exceptions


# instrument and radius of search in arsec.
field_radius = {"uvit": 1200, "sxt": 1500, "czti": 1680, "laxpc": 1680}

# instrument and required window size.
field_radius_im = {"uvit": 800, "sxt": 1500, "czti": 1120, "laxpc": 1120}

# TD1 related.
flux_norm = 2e-13


# To convert B-V to Spectral Type (reference: http://www.stsci.edu/~inr/intrins.html)
def spectype(f):
    if f >= 1.511:
        return "M4"
    elif f >= 1.486:
        return "M3"
    elif f >= 1.421:
        return "M1"
    elif f >= 1.351:
        return "M0"
    elif f >= 1.241:
        return "K7"
    elif f >= 1.076:
        return "K5"
    elif f >= 0.976:
        return "K4"
    elif f >= 0.936:
        return "K3"
    elif f >= 0.891:
        return "K2"
    elif f >= 0.836:
        return "K1"
    elif f >= 0.776:
        return "K0"
    elif f >= 0.711:
        return "G8"
    elif f >= 0.666:
        return "G5"
    elif f >= 0.641:
        return "G3"
    elif f >= 0.616:
        return "G2"
    elif f >= 0.566:
        return "G0"
    elif f >= 0.491:
        return "F8"
    elif f >= 0.401:
        return "F5"
    elif f >= 0.346:
        return "F2"
    elif f >= 0.331:
        return "F1"
    elif f >= 0.311:
        return "F0"
    elif f >= 0.286:
        return "A9"
    elif f >= 0.236:
        return "A8"
    elif f >= 0.186:
        return "A7"
    elif f >= 0.161:
        return "A6"
    elif f >= 0.136:
        return "A5"
    elif f >= 0.101:
        return "A4"
    elif f >= 0.066:
        return "A3"
    elif f >= 0.036:
        return "A2"
    elif f >= 0.006:
        return "A1"
    elif f >= -0.039:
        return "A0"
    elif f >= -0.089:
        return "B9"
    elif f >= -0.119:
        return "B8"
    elif f >= -0.134:
        return "B7"
    elif f >= -0.149:
        return "B6"
    elif f >= -0.169:
        return "B5"
    elif f >= -0.189:
        return "B4"
    elif f >= -0.219:
        return "B3"
    elif f >= -0.249:
        return "B2"
    elif f >= -0.279:
        return "B1"
    elif f < -0.279:
        return "B0"
    else:
        print("out of bounds")


def observe_VIS(instrument, RA, DEC):
    # This code is going to throw some warnings.
    warnings.filterwarnings("ignore")

    # proximity parameter.
    proximity = 10.0

    # Checking if the UVIT server is live.
    uvit_session = Session()

    # To run beloved BSWT.
    RADEC = RA.replace(":", " ") + ", " + DEC.replace(":", " ")
    bswt_data = {
        "coord_type": "eq",
        "coords": RADEC,
        "source": "",
        "prinst": instrument,
    }

    try_value = 0
    got_response = False
    while got_response is False:
        try:
            bswt_html = uvit_session.post(
                url="https://uvit.iiap.res.in/cgi-bin/bswt.pl",
                data=bswt_data,
                verify=False,
            )
            got_response = True
        except (exceptions.RequestException, AttributeError):
            try_value = try_value + 1
            if try_value == 50:
                print(
                    """\nCanUVIT attempted access to the BSWT online tool 50 times without success. Please check:

    1. Your internet connection.
    2. Accessibility of uvit.iiap.res.in."""
                )
                return None

    bswt_soup = BeautifulSoup(bswt_html.text, "html.parser")
    text_output_link = bswt_soup.find("a", href=True)
    text_output = uvit_session.get(
        "https://uvit.iiap.res.in" + text_output_link["href"]
    )
    bswt_data = np.genfromtxt(
        BytesIO(bytes(text_output.text, "UTF-8")), skip_header=10, invalid_raise=False
    )

    bswt_data = np.unique(bswt_data[:, 2:-2:2], axis=0)

    # To get the magnitude and color
    if len(bswt_data.shape) == 2:
        c1l, c2l, c3l, c4l = zip(*bswt_data)
        sortedl = sorted(zip(map(float, c3l), c4l, c1l, c2l))
        cm, ct, ra_deg, dec_deg = zip(*sortedl)
    else:
        c1l, c2l, c3l, c4l = bswt_data
        sortedl = [float(c3l), c4l, c1l, c2l]
        sortedl = [[i] for i in sortedl]
        cm, ct, ra_deg, dec_deg = sortedl

    cm = cm[0:7]
    ct = ct[0:7]
    ra_deg = ra_deg[0:7]
    dec_deg = dec_deg[0:7]
    fct = map(float, list(ct))
    spty = list(map(spectype, fct))

    # iteratively inputing parameters to etc.
    vs3list = []
    vs2list = []
    vs1list = []
    nd1list = []
    bk7list = []
    for i in range(len(cm)):
        # input parameters
        src_mag = str(cm[i])
        sptype1 = spty[i][0]
        sptype2 = spty[i][-1]

        # form encoded data
        etcdata = {
            "src_type": "star",
            "sptype1": sptype1,
            "sptype2": sptype2,
            "sptype3": "V",
            "bbodytemp": "6000.0",
            "galaxyclass": "sc",
            "agnclass": "seyfert2",
            "plaw_index": "-1.0",
            "fluxval": "2.0",
            "flatspec_unit": "cgs",
            "redshift": "0.00",
            "ftype": "usemag",
            "src_mag": src_mag,
            "mag_band": "v",
            "coords": "11 00 00.00, -16 00 00.0",
            "ctype": "equatorial",
            "ra": "0",
            "dec": "0",
            "rv": "3.1",
            "ebv": "0.0",
            "nh": "1.00",
            "distance": "0.45",
            "av": "1.0",
            "ext_mode": "rvebv",
            "dc": "25",
            "calc": "et",
            "snr": "5.0",
            "et": "1800",
        }

        # magic!
        session = Session()
        try_value = 0
        got_response = False
        while got_response is False:
            try:
                response = session.post(
                    url="https://uvit.iiap.res.in/cgi-bin/etc.pl",
                    data=etcdata,
                    verify=False,
                )

                got_response = True
            except (exceptions.RequestException, AttributeError):
                try_value = try_value + 1
                if try_value == 50:
                    print(
                        """\nCanUVIT attempted access to the ETC online tool 50 times without success. Please check:

    1. Your internet connection.
    2. Accessibility of uvit.iiap.res.in."""
                    )
                    return None

        # retrieving required values using soup
        soup = BeautifulSoup(response.text, "html.parser")
        ss = soup.find("table", {"id": "aux"})
        rows = ss.find_all("tr")
        for row in rows:
            if row.find("td").get_text() == "VIS 3":
                vis3 = row.find("td").next_element.next_element.get_text()

                if (
                    len(vis3) >= 10
                ):  # A work-around (example: turns "1.40 x 10+04" to "14000.0")
                    vis3l = vis3.split()
                    vis3 = float(vis3l[0]) * (10 ** float(vis3l[2][-1]))

                vis3 = float(vis3)
                vs3list.append(vis3)

            if row.find("td").get_text() == "VIS 2":
                vis2 = row.find("td").next_element.next_element.get_text()

                if len(vis2) >= 10:
                    vis2l = vis2.split()
                    vis2 = float(vis2l[0]) * (10 ** float(vis2l[2][-1]))

                vis2 = float(vis2)
                vs2list.append(vis2)

            if row.find("td").get_text() == "VIS 1":
                vis1 = row.find("td").next_element.next_element.get_text()

                if len(vis1) >= 10:
                    vis1l = vis1.split()
                    vis1 = float(vis1l[0]) * (10 ** float(vis1l[2][-1]))

                vis1 = float(vis1)
                vs1list.append(vis1)

            if row.find("td").get_text() == "VIS ND1":
                nd1 = row.find("td").next_element.next_element.get_text()

                if len(nd1) >= 10:
                    ndl = nd1.split()
                    nd1 = float(ndl[0]) * (10 ** float(ndl[2][-1]))

                nd1 = float(nd1)
                nd1list.append(nd1)

            if row.find("td").get_text() == "VIS BK-7":
                bk7 = row.find("td").next_element.next_element.get_text()

                if len(bk7) >= 10:
                    bkl = bk7.split()
                    bk7 = float(bkl[0]) * (10 ** float(bkl[2][-1]))

                bk7 = float(bk7)
                bk7list.append(bk7)

    # To convert ra_deg and dec_deg to ra_hms and dec_dms.
    coord = SkyCoord(list(zip(ra_deg, dec_deg)), frame="icrs", unit="deg")
    RAhms_DECdms = coord.to_string("hmsdms", sep=":", precision=1)
    ra_hms, dec_dms = zip(*[hmdm.split(" ") for hmdm in RAhms_DECdms])

    # Conversion to numpy arrays.
    ra_deg = np.array(ra_deg)
    dec_deg = np.array(dec_deg)
    ra_hms = np.array(ra_hms)
    dec_dms = np.array(dec_dms)
    cm = np.array(cm)
    ct = np.array(ct)
    spty = np.array(spty)
    vs3list = np.array(vs3list)
    vs2list = np.array(vs2list)
    vs1list = np.array(vs1list)
    nd1list = np.array(nd1list)
    bk7list = np.array(bk7list)

    vs3list = np.round(vs3list, 1)
    vs2list = np.round(vs2list, 1)
    vs1list = np.round(vs1list, 1)
    nd1list = np.round(nd1list, 1)
    bk7list = np.round(bk7list, 1)

    sl_no = np.arange(len(ra_hms)) + 1

    vis_res = Table(
        [
            sl_no,
            ra_hms,
            dec_dms,
            cm,
            ct,
            spty,
            vs3list,
            vs2list,
            vs1list,
            nd1list,
            bk7list,
        ],
        names=(
            "sl_no",
            "ra_hms",
            "dec_dms",
            "mag",
            "B-V",
            "SpecType",
            "VIS3",
            "VIS2",
            "VIS1",
            "ND1",
            "BK7",
        ),
        meta={"name": "VIS counts"},
    )

    vis_res["mag"].format = ".1f"
    vis_res["B-V"].format = ".1f"
    vis_res["VIS3"].format = ".1f"
    vis_res["VIS2"].format = ".1f"
    vis_res["VIS1"].format = ".1f"
    vis_res["ND1"].format = ".1f"
    vis_res["BK7"].format = ".1f"

    # To select safe filters.
    filter_dict = {0: "VIS3", 1: "VIS2", 2: "VIS1", 3: "ND1", 4: "BK7"}
    i = 0
    safe_filters = []
    for Filter in [vs3list, vs2list, vs1list, nd1list, bk7list]:
        if sum(Filter > 4800) == 0:
            safe_filters.append(filter_dict[i])
        i = i + 1

    # To check if objects within 10 arcsec are present.
    sep = [
        ca.separation(cb).arcsecond
        for ca in coord
        for cb in coord
        if ca.ra.value != cb.ra.value
    ]

    proximity_check = np.array(sep) < proximity
    too_close = sum(proximity_check) / 2

    print("\n### VIS\n")
    vis_res.pprint_all()

    print("\nSafe filters in VIS: {}".format(safe_filters))

    if too_close > 0:
        print(
            "\nWARNING! There exists {} pair of bright stars which are closer than\
              \n{} arcseconds!".format(
                too_close, proximity
            )
        )


# Function to find seperation in celestial coordinates.
def cel_separation(RA_deg, DEC_deg, pointing_coo):
    coo = SkyCoord(RA_deg, DEC_deg, frame="icrs", unit="deg")
    return coo.separation(pointing_coo)


# Function to convert ra_deg and dec_deg to ra_hms and dec_dms.
def deg_to_hms(ra_deg, dec_deg):
    fuv_coord = SkyCoord(np.array([ra_deg, dec_deg]).T, frame="icrs", unit="deg")

    ra_hms_dec_dms = fuv_coord.to_string("hmsdms", sep=":", precision=1)
    ra_hms, dec_dms = list(zip(*[hmdm.split(" ") for hmdm in ra_hms_dec_dms]))
    sl_no = np.arange(len(ra_hms)) + 1
    xy_tab = Table([sl_no, ra_hms, dec_dms], names=("sl_no", "ra_hms", "dec_dms"))

    return xy_tab


# Functions to convert GALEX CPS to AB magnitude.
def magfuv(FUV_CPS):
    FUV_ZPMAG = 18.82
    return -2.5 * np.log10(np.array(FUV_CPS)) + FUV_ZPMAG


def magnuv(NUV_CPS):
    NUV_ZPMAG = 20.08
    return -2.5 * np.log10(np.array(NUV_CPS)) + NUV_ZPMAG


# Functions to convert GALEX magnitude to UVIT count rates.
def countfuv(MAG):
    caf2 = 1.0
    baf2 = 0.85
    sapphire = 0.63
    silica = 0.22
    MAG1 = 18.22

    if 10.511 <= MAG <= 15.0:
        MAG_c = 5.371 + (20.0 * MAG - 210.2) ** 0.5
    elif MAG < 10.511:
        MAG_c = 5.371
    else:
        MAG_c = MAG

    cr1 = caf2 * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr2 = baf2 * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr3 = sapphire * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr4 = silica * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    return MAG, MAG_c, cr1, cr2, cr3, cr4


countfuv = np.vectorize(countfuv)


def countfuv_abs(MAG):  # for cases where FUV is absent.
    caf2 = 1.0
    baf2 = 0.85
    sapphire = 0.63
    silica = 0.22
    MAG1 = 18.22

    if 9.323 <= MAG <= 15.0:
        MAG_c = 2.634 + (26.316 * MAG - 245.329) ** 0.5
    elif MAG < 9.323:
        MAG_c = 2.634
    else:
        MAG_c = MAG

    MAG_c = MAG_c - 1.65
    cr1 = caf2 * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr2 = baf2 * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr3 = sapphire * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr4 = silica * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    return MAG, MAG_c, cr1, cr2, cr3, cr4


countfuv_abs = np.vectorize(countfuv_abs)


def countnuv(MAG):
    silica = 1.0
    b4 = 0.22
    b13 = 0.27
    b15 = 0.074
    n2 = 0.055
    MAG1 = 20.0

    if 9.323 <= MAG <= 15.0:
        MAG_c = 2.634 + (26.316 * MAG - 245.329) ** 0.5
    elif MAG < 9.323:
        MAG_c = 2.634
    else:
        MAG_c = MAG

    cr1 = silica * 10.0 ** ((MAG1 - MAG_c) * 0.4)
    cr2 = b4 * 10 ** ((MAG1 - MAG_c) * 0.4)
    cr3 = b13 * 10 ** ((MAG1 - MAG_c) * 0.4)
    cr4 = b15 * 10 ** ((MAG1 - MAG_c) * 0.4)
    cr5 = n2 * 10 ** ((MAG1 - MAG_c) * 0.4)
    return MAG, MAG_c, cr1, cr2, cr3, cr4, cr5


countnuv = np.vectorize(countnuv)


# Function to detect blobs, estimate fluxes, and do some more.
def im_flux(int_map, pointing_coo, instrument, m_ra, m_dec):
    # To read data.
    try:
        hdu = fits.open(int_map, cache=False)
    except IOError:
        print("Incomplete FITS file. Please try again.")

    # To convert RA & DEC to pixel coordinates.
    w = WCS(hdu[0].header)
    cor = w.all_world2pix(pointing_coo.ra.degree, pointing_coo.dec.degree, 1)
    try:
        selcen = [int(np.round(x)) for x in cor]
    except ValueError:
        print("The provided RA DEC values fell outside the GALEX image.")

    fitsf = hdu[0].data

    # The pixel scale of GALEX taken here is 1.5 arcsec/pixel.
    # To select a rectangular region centered on the provided positions.
    xfi = int(selcen[0]) - field_radius_im[instrument]
    yfi = int(selcen[1]) - field_radius_im[instrument]
    xse = int(selcen[0]) + field_radius_im[instrument]
    yse = int(selcen[1]) + field_radius_im[instrument]
    xfi = max(0, xfi)
    yfi = max(0, yfi)
    fitsf = fitsf[yfi:yse, xfi:xse]

    # To put a circular mask on the image.
    yo, xo = np.ogrid[
        -field_radius_im[instrument] : field_radius_im[instrument],
        -field_radius_im[instrument] : field_radius_im[instrument],
    ]

    mask = xo * xo + yo * yo > field_radius_im[instrument] ** 2

    # To see if mask needs to be reduced in size.
    mask_xshape, mask_yshape = mask.shape
    fitsf_xshape, fitsf_yshape = fitsf.shape
    x_mask_start = mask_xshape - fitsf_xshape
    y_mask_start = mask_yshape - fitsf_yshape
    mask = mask[x_mask_start:, y_mask_start:]
    fitsf[mask] = 0

    # To mark positions of bright objects
    det_pos = w.all_world2pix(m_ra, m_dec, 1)
    x_det = det_pos[0] - xfi
    y_det = det_pos[1] - yfi

    # Using a 7x7 box to estimate the flux of detected objects.
    # the units are counts/sec/pixel.
    fluxes = []
    for blob in zip(y_det, x_det):
        ydet, xdet = blob
        xfirst = int(xdet) - 3
        yfirst = int(ydet) - 3
        xsecond = int(xdet) + 4
        ysecond = int(ydet) + 4
        fluxes.append(fitsf[yfirst:ysecond, xfirst:xsecond].sum())

    # To plot the image.
    plt.imshow(fitsf, cmap="gray", norm=LogNorm(), interpolation="none")

    plt.title("Detected bright sources marked")
    plt.gca().invert_yaxis()
    plt.xticks([])
    plt.yticks([])

    # To plot detected sources.
    plt.scatter(x_det, y_det, color="r", marker="o", alpha=0.2)

    # To annotate positions.
    anno = np.arange(len(x_det)) + 1
    for q, txt in enumerate(anno):
        plt.annotate(txt, (x_det[q], y_det[q]))

    # To save the image.
    figure_name = int_map.split("/")[-1].replace(".fits.gz", ".png")
    plt.savefig(figure_name, format="png", bbox_inches="tight", dpi=300)

    plt.clf()
    return fluxes


# Functions to TD1 flux to UVIT count rates.
def td1_countnuv(flux):
    silica = 955.0
    b4 = 218.5
    b13 = 275.8
    b15 = 59.6
    n2 = 50.6
    flux_ratio = flux / flux_norm
    cr1 = silica * flux_ratio
    cr2 = b4 * flux_ratio
    cr3 = b13 * flux_ratio
    cr4 = b15 * flux_ratio
    cr5 = n2 * flux_ratio
    return flux, cr1, cr2, cr3, cr4, cr5


def td1_countfuv(flux):
    caf2 = 74.5
    baf2 = 60.0
    sapphire = 50.0
    silica = 17.3
    flux_ratio = flux / flux_norm
    cr1 = caf2 * flux_ratio
    cr2 = baf2 * flux_ratio
    cr3 = sapphire * flux_ratio
    cr4 = silica * flux_ratio
    return flux, cr1, cr2, cr3, cr4


# Function to do all the work on TD1_catalogue.
def td1_estimate(pointing_coo, instrument):
    td1_catalogue = "https://github.com/prajwel/canuvit/blob/main/src/canuvit/td1_catalogue.fits?raw=true"
    td1_hdu = fits.open(td1_catalogue, cache=False)
    alpha = td1_hdu[1].data["ra"]
    delta = td1_hdu[1].data["dec"]
    nuv_flux = td1_hdu[1].data["flux_2365_a"]
    fuv_flux = td1_hdu[1].data["flux_1565_a"]

    # NUV
    refined_set = [
        (al, de, nf)
        for al, de, nf in zip(alpha, delta, nuv_flux)
        if (pointing_coo.ra.value - 5) <= al <= (pointing_coo.ra.value + 5)
        and (pointing_coo.dec.value - 5) <= de <= (pointing_coo.dec.value + 5)
    ]

    nalpha, ndelta, nuv_flux = list(zip(*refined_set))

    confined_set = [
        nf
        for al, de, nf in zip(nalpha, ndelta, nuv_flux)
        if cel_separation(al, de, pointing_coo) <= field_radius[instrument] * u.arcsec
    ]

    # If list is empty, normal value need to be taken.
    if len(confined_set) == 0:
        confined_set.append(flux_norm)

    nd = sorted(confined_set)[-1]
    flux, ta, tb, tc, td, te = td1_countnuv(nd)
    nuv_res = Table(
        [[flux], [ta], [tb], [tc], [td], [te]],
        names=("flux_2365_a", "silica", "b4", "b13", "b15", "n2"),
        meta={"name": "NUV counts"},
    )

    nuv_res["silica"].format = ".1f"
    nuv_res["b4"].format = ".1f"
    nuv_res["b13"].format = ".1f"
    nuv_res["b15"].format = ".1f"
    nuv_res["n2"].format = ".1f"

    print("\n### NUV\n")
    nuv_res.pprint_all()

    # To select NUV safe filters.
    nuv_filter_dict = {
        0: "Silica",
        1: "NUV-B4",
        2: "NUV-B13",
        3: "NUV-B15",
        4: "NUV-N2",
    }
    i = 0
    nuv_safe = []
    for Filter in list(zip(*nuv_res["silica", "b4", "b13", "b15", "n2"])):
        if sum(np.array(Filter) > 1500) == 0:
            nuv_safe.append(nuv_filter_dict[i])
        if i == 0:
            if sum(np.array(Filter) > 1133) == 0:
                nuv_safe.append("NUV-grating")
        i = i + 1

    nuv_declaration = "Safe filters in NUV: {}".format(nuv_safe)
    print("\n{}".format(nuv_declaration))

    # FUV
    refined_set = [
        (al, de, ff)
        for al, de, ff in zip(alpha, delta, fuv_flux)
        if (pointing_coo.ra.value - 5) <= al <= (pointing_coo.ra.value + 5)
        and (pointing_coo.dec.value - 5) <= de <= (pointing_coo.dec.value + 5)
    ]

    nalpha, ndelta, fuv_flux = list(zip(*refined_set))

    confined_set = [
        ff
        for al, de, ff in zip(nalpha, ndelta, fuv_flux)
        if cel_separation(al, de, pointing_coo) <= field_radius[instrument] * u.arcsec
    ]

    # If list is empty, normal value need to be taken.
    if len(confined_set) == 0:
        confined_set.append(flux_norm)

    fd = sorted(confined_set)[-1]
    flux, ta, tb, tc, td = td1_countfuv(fd)
    fuv_res = Table(
        [[flux], [ta], [tb], [tc], [td]],
        names=("flux_1565_a", "caf2", "baf2", "sapphire", "silica"),
        meta={"name": "NUV counts"},
    )

    fuv_res["caf2"].format = ".1f"
    fuv_res["baf2"].format = ".1f"
    fuv_res["sapphire"].format = ".1f"
    fuv_res["silica"].format = ".1f"

    print("\n### FUV\n")
    fuv_res.pprint_all()

    # To select FUV safe filters.
    fuv_filter_dict = {0: "CaF2", 1: "BaF2", 2: "Sapphire", 3: "Silica"}
    j = 0
    fuv_safe = []
    for Filter in list(zip(*fuv_res["caf2", "baf2", "sapphire", "silica"])):
        if sum(np.array(Filter) > 1500) == 0:
            fuv_safe.append(fuv_filter_dict[j])
        if j == 0:
            if sum(np.array(Filter) > 892) == 0:
                fuv_safe.append("FUV-grating")
        j = j + 1

    fuv_declaration = "Safe filters in FUV: {}".format(fuv_safe)
    print("\n{}".format(fuv_declaration))


# Function to format NUV data.
def format_nuv(nuv_counts):
    ntab = Table(
        nuv_counts,
        names=("Mag", "Mag_corrected", "silica", "b4", "b13", "b15", "n2"),
        meta={"name": "NUV counts"},
    )

    ntab["Mag"].format = ".1f"
    ntab["Mag_corrected"].format = ".1f"
    ntab["silica"].format = ".1f"
    ntab["b4"].format = ".1f"
    ntab["b13"].format = ".1f"
    ntab["b15"].format = ".1f"
    ntab["n2"].format = ".1f"
    return ntab


# Function to format FUV data.
def format_fuv(fuv_counts):
    ftab = Table(
        fuv_counts,
        names=("Mag", "Mag_corrected", "caf2", "baf2", "sapphire", "silica"),
        meta={"name": "FUV counts"},
    )

    ftab["Mag"].format = ".1f"
    ftab["Mag_corrected"].format = ".1f"
    ftab["caf2"].format = ".1f"
    ftab["baf2"].format = ".1f"
    ftab["sapphire"].format = ".1f"
    ftab["silica"].format = ".1f"
    return ftab


def observe_UV(instrument, RA, DEC):
    pointing_coo = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg))
    instrument = instrument.lower()

    # Check if source Galactic latitude is between -30 to 30.
    gal_lat = pointing_coo.galactic.b.value
    gal_plane = "no"
    if -30.0 <= gal_lat <= 30.0:
        gal_plane = "yes"

    # Get on with MAST website queries.
    obs_table = Observations.query_region(
        str(pointing_coo.ra.degree) + " " + str(pointing_coo.dec.degree),
        radius=".001 deg",
    )

    obs_table = obs_table[obs_table["obs_collection"] == "GALEX"]
    obs_table = obs_table[obs_table["dataproduct_type"] == "image"]

    # To make sense of the mess that is MAST.
    if len(obs_table) == 0:
        no_galex_tiles = "0 Galex tiles found. Galex observations around \
                          \nthe given target is not available. Using TD1\
                          \ncatalogue to estimate UVIT count rates."
        print("\n{}".format(no_galex_tiles))
        if gal_plane == "yes":
            gal_plane_warning = "The galactic latitude is between -30 to 30. \
                                \nYour field cannot be checked using TD1 catalogue!"
            print("\n{}".format(gal_plane_warning))
            return
        else:
            # To read the TD1 catalogue.
            td1_estimate(pointing_coo, instrument)
            return

    if len(obs_table) > 1:
        obs_table_coo = SkyCoord(obs_table["s_ra"], obs_table["s_dec"], unit=u.deg)
        obs_table_offsets = pointing_coo.separation(obs_table_coo)
        sort_indices = np.argsort(obs_table_offsets.degree)
        obs_table = obs_table[sort_indices]

    data_products = Observations.get_product_list(obs_table[0])
    data_products = data_products[
        data_products["productGroupDescription"] == "Minimum Recommended Products"
    ]

    catalogue_link = []
    image_links = []
    base_url = "https://mast.stsci.edu/api/v0.1/Download/file?uri="
    for somel in data_products["dataURI"]:
        try:
            if somel[-12:] == "mcat.fits.gz":
                catalogue_link.append(base_url + somel)
            elif somel[-11:] == "int.fits.gz":
                image_links.append(base_url + somel)
        except Exception:
            pass

    if len(catalogue_link) != 0:
        try:
            cat_hdu = fits.open(catalogue_link[0], cache=False)
        except IOError:
            incomplete_fits = "Incomplete FITS file. Please try again."
            print("\n{}".format(incomplete_fits))
    else:
        no_catalogue = "Could not find the catalogue for this region."
        print("\n{}".format(no_catalogue))

    alpha = cat_hdu[1].data["alpha_j2000_merged"]
    delta = cat_hdu[1].data["delta_j2000_merged"]

    # NUV
    nuv_mag = cat_hdu[1].data["nuv_mag"]
    nuv_fwhm = cat_hdu[1].data["nuv_fwhm_world"]
    refined_set = [
        (al, de, nm, nf)
        for al, de, nm, nf in zip(alpha, delta, nuv_mag, nuv_fwhm)
        if int(nm) != -999 and nm <= 22.0
    ]

    nalpha, ndelta, nuv_mag, nuv_fwhm = list(zip(*refined_set))

    confined_set = [
        (nm, al, de, nf)
        for al, de, nm, nf in zip(nalpha, ndelta, nuv_mag, nuv_fwhm)
        if cel_separation(al, de, pointing_coo) <= field_radius[instrument] * u.arcsec
    ]

    nd = np.array(sorted(confined_set))[0:5]
    cat_nuv_counts = countnuv(nd[:, 0])
    cat_nuv_res = format_nuv(cat_nuv_counts)

    # To convert ra_deg and dec_deg to ra_hms and dec_dms.
    xy_tab = deg_to_hms(nd[:, 1], nd[:, 2])
    cat_nuv_res = hstack([xy_tab, cat_nuv_res])

    balance = cat_nuv_res["ra_hms", "dec_dms", "Mag"]
    balance.rename_column("Mag", "CAT_Mag")
    balance["fwhm"] = nd[:, 3]
    balance["fwhm"].format = "4.4f"

    # FUV
    fuv_mag = cat_hdu[1].data["fuv_mag"]

    fuv_absent = "no"
    if len(np.unique(fuv_mag)) == 1:  # when FUV data is absent.
        fd = nd
        warning = "\nFUV observations seem to be absent! Using M_fuv = M_nuv - 1.65."
        print(warning)
        fuv_absent = "yes"
    else:
        refined_set = [
            (al, de, fm)
            for al, de, fm in zip(alpha, delta, fuv_mag)
            if int(fm) != -999 and fm <= 22.0
        ]

        falpha, fdelta, fuv_mag = list(zip(*refined_set))

        confined_set = [
            (fm, al, de)
            for al, de, fm in zip(falpha, fdelta, fuv_mag)
            if cel_separation(al, de, pointing_coo)
            <= field_radius[instrument] * u.arcsec
        ]

        fd = np.array(sorted(confined_set))[0:5]

    if fuv_absent == "no":
        cat_fuv_counts = countfuv(fd[:, 0])
    else:
        cat_fuv_counts = countfuv_abs(fd[:, 0])

    cat_fuv_res = format_fuv(cat_fuv_counts)

    # To convert ra_deg and dec_deg to ra_hms and dec_dms.
    xy_tab = deg_to_hms(fd[:, 1], fd[:, 2])
    cat_fuv_res = hstack([xy_tab, cat_fuv_res])

    if len(image_links) != 0:
        for image_link in image_links:
            detector = str(image_link[-14:-13]).lower()
            if detector == "f":
                fuv_intmap = image_link
            elif detector == "n":
                nuv_intmap = image_link
            else:
                filter_confusion = "Cannot determine filter! Exiting."
                print("\n{}".format(filter_confusion))

    # NUV
    if "nuv_intmap" in locals():
        m_ra = nd[:, 1]
        m_dec = nd[:, 2]
        n_fluxes = im_flux(nuv_intmap, pointing_coo, instrument, m_ra, m_dec)
        n_mags = magnuv(n_fluxes)
        n_counts = countnuv(n_mags)
        im_nuv_res = format_nuv(n_counts)
        xy_tab = deg_to_hms(m_ra, m_dec)
        im_nuv_res = hstack([xy_tab, im_nuv_res])
        balance["IM_Mag"] = im_nuv_res["Mag"]

    # FUV
    if "fuv_intmap" not in locals() and "nuv_intmap" in locals():
        f_counts = countfuv_abs(n_mags)
        im_fuv_res = format_fuv(f_counts)
        im_fuv_res = hstack([xy_tab, im_fuv_res])

    elif "fuv_intmap" in locals():
        m_ra = fd[:, 1]
        m_dec = fd[:, 2]
        f_fluxes = im_flux(fuv_intmap, pointing_coo, instrument, m_ra, m_dec)
        f_mags = magfuv(f_fluxes)
        f_counts = countfuv(f_mags)
        im_fuv_res = format_fuv(f_counts)
        xy_tab = deg_to_hms(m_ra, m_dec)
        im_fuv_res = hstack([xy_tab, im_fuv_res])

    # To decide between catalogue or image.
    balance["diff"] = balance["IM_Mag"] - balance["CAT_Mag"]

    if sum(balance["diff"] > 1.2) == 0 and sum(balance["fwhm"] > 0.0043) == 0:
        nuv_res = cat_nuv_res
        fuv_res = cat_fuv_res
    else:
        nuv_res = im_nuv_res
        fuv_res = im_fuv_res

    # To select NUV safe filters.
    print("\n### NUV\n")
    nuv_res.pprint_all()
    nuv_filter_dict = {
        0: "Silica",
        1: "NUV-B4",
        2: "NUV-B13",
        3: "NUV-B15",
        4: "NUV-N2",
    }
    i = 0
    nuv_safe = []
    for Filter in zip(*nuv_res["silica", "b4", "b13", "b15", "n2"]):
        if sum(np.array(Filter) > 1500) == 0:
            nuv_safe.append(nuv_filter_dict[i])
        if i == 0:
            if sum(np.array(Filter) > 1133) == 0:
                nuv_safe.append("NUV-grating")
        i = i + 1

    nuv_declaration = "Safe filters in NUV: {}".format(nuv_safe)
    print("\n{}".format(nuv_declaration))

    # To select FUV safe filters.
    print("\n### FUV\n")
    fuv_res.pprint_all()
    fuv_filter_dict = {0: "CaF2", 1: "BaF2", 2: "Sapphire", 3: "Silica"}
    j = 0
    fuv_safe = []
    for Filter in zip(*fuv_res["caf2", "baf2", "sapphire", "silica"]):
        if sum(np.array(Filter) > 1500) == 0:
            fuv_safe.append(fuv_filter_dict[j])
        if j == 0:
            if sum(np.array(Filter) > 892) == 0:
                fuv_safe.append("FUV-grating")
        j = j + 1

    fuv_declaration = "Safe filters in FUV: {}".format(fuv_safe)
    print("\n{}".format(fuv_declaration))


def observe(instrument, RA, DEC):
    print("\n=========================================================")
    print("Payload: {}, Coordinates: {}, {}".format(instrument, RA, DEC))
    print("=========================================================")

    observe_VIS(instrument, RA, DEC)
    observe_UV(instrument, RA, DEC)
