#!/usr/bin/env python3


import sys
import warnings
import numpy as np
from io import BytesIO

from bs4 import BeautifulSoup
#from astropy import units as u
from requests import Session, exceptions
from astropy.coordinates import SkyCoord


# To convert B-V to Spectral Type (reference: http://www.stsci.edu/~inr/intrins.html) 
def spectype(f):
    if f >= 1.511:
        return "M4"
    elif f >= 1.486:
        return  "M3"
    elif f >= 1.421:
        return  "M1"
    elif f >= 1.351: 
        return  "M0"
    elif f >= 1.241: 
        return  "K7"
    elif f >= 1.076: 
        return  "K5"
    elif f >= 0.976: 
        return  "K4"
    elif f >= 0.936: 
        return  "K3"
    elif f >= 0.891: 
        return  "K2"
    elif f >= 0.836: 
        return  "K1"
    elif f >= 0.776: 
        return  "K0"
    elif f >= 0.711: 
        return  "G8"
    elif f >= 0.666: 
        return  "G5"
    elif f >= 0.641: 
        return  "G3"
    elif f >= 0.616: 
        return  "G2"
    elif f >= 0.566: 
        return  "G0"
    elif f >= 0.491: 
        return  "F8"
    elif f >= 0.401: 
        return  "F5"
    elif f >= 0.346: 
        return  "F2"
    elif f >= 0.331: 
        return  "F1"
    elif f >= 0.311: 
        return  "F0"
    elif f >= 0.286: 
        return  "A9"
    elif f >= 0.236: 
        return  "A8"
    elif f >= 0.186: 
        return  "A7"
    elif f >= 0.161: 
        return  "A6"
    elif f >= 0.136: 
        return  "A5"
    elif f >= 0.101: 
        return  "A4"
    elif f >= 0.066: 
        return  "A3"
    elif f >= 0.036: 
        return  "A2"
    elif f >= 0.006: 
        return  "A1"
    elif f >= -0.039:
        return  "A0"
    elif f >= -0.089:
        return  "B9"
    elif f >= -0.119:
        return  "B8"
    elif f >= -0.134:
        return  "B7"
    elif f >= -0.149:
        return  "B6"
    elif f >= -0.169:
        return  "B5"
    elif f >= -0.189:
        return  "B4"
    elif f >= -0.219:
        return  "B3"
    elif f >= -0.249:
        return  "B2"
    elif f >= -0.279:
        return  "B1"
    elif f < -0.279:
        return  "B0"
    else: 
        print("out of bounds")
        
def theia(instrument, RA_user, DEC_user):

    # This code is going to throw some warnings.
    warnings.filterwarnings("ignore")

    # proximity parameter.
    proximity = 10.
        
    # Checking if the UVIT server is live.
    uvit_session = Session()
        
    # To run beloved BSWT.
    RADEC = RA_user.replace(':', ' ') + ', ' + DEC_user.replace(':', ' ')
    bswt_data = {'coord_type': 'eq',
                'coords': RADEC,
                'source': '',
                'prinst': instrument}

    try_value = 0
    got_response = False
    while got_response == False:
        try:
            bswt_html = uvit_session.post(url = 'https://uvit.iiap.res.in/cgi-bin/bswt.pl',
                                          data = bswt_data, verify = False)   
            got_response = True
        except (exceptions.RequestException, AttributeError):
            try_value = try_value + 1
            if try_value == 50:
                print('Check your network, Theia tried 50 times without a luck')
                return None

    bswt_soup = BeautifulSoup(bswt_html.text, 'html.parser')
    text_output_link = bswt_soup.find('a', href = True)
    text_output = uvit_session.get("https://uvit.iiap.res.in" + text_output_link['href'])
    bswt_data = np.genfromtxt(BytesIO(bytes(text_output.text, 'UTF-8')), 
                              skip_header = 10, 
                              invalid_raise = False)

    bswt_data = np.unique(bswt_data[:, 2: -2: 2], axis= 0)

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

    cm = cm[0: 7] 
    ct = ct[0: 7] 
    ra_deg = ra_deg[0: 7]
    dec_deg = dec_deg[0: 7]
    fct = map(float, list(ct))
    spty = list(map(spectype, fct))

    #iteratively inputing parameters to etc. 
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

        #form encoded data
        etcdata = {'src_type': 'star', 'sptype1': sptype1, 'sptype2': sptype2, 
                   'sptype3': 'V', 'bbodytemp': '6000.0', 'galaxyclass': 'sc', 
                   'agnclass': 'seyfert2', 'plaw_index': '-1.0', 'fluxval': '2.0', 
                   'flatspec_unit': 'cgs', 'redshift': '0.00', 'ftype': 'usemag', 
                   'src_mag': src_mag, 'mag_band': 'v', 'coords': '11 00 00.00, -16 00 00.0', 
                   'ctype': 'equatorial', 'ra': '0', 'dec': '0', 'rv': '3.1', 
                   'ebv': '0.0', 'nh': '1.00', 'distance': '0.45', 'av': '1.0', 
                   'ext_mode': 'rvebv', 'dc': '25', 'calc': 'et', 'snr': '5.0', 'et': '1800'}
        
        #magic!
        session = Session()
        try_value = 0
        got_response = False
        while got_response == False:
            try:
                response = session.post(url = "https://uvit.iiap.res.in/cgi-bin/etc.pl", 
                                        data = etcdata, 
                                        verify = False)   

                got_response = True
            except (exceptions.RequestException, AttributeError):
                try_value = try_value + 1
                if try_value == 50:
                    print('Check your network, Theia tried 50 times without a luck')
                    return None

        #retrieving required values using soup
        soup = BeautifulSoup(response.text, 'html.parser')
        ss = soup.find("table", { "id" : "aux" }) 
        rows = ss.find_all('tr')    
        for row in rows:
        
            if row.find('td').get_text() == 'VIS 3':
                vis3 = row.find('td').next_element.next_element.get_text()
        
                if len(vis3) >= 10: # A work-around (example: turns "1.40 x 10+04" to "14000.0")  
                    vis3l = vis3.split()
                    vis3 = float(vis3l[0]) * (10 ** float(vis3l[2][-1]))
        
                vis3 = float(vis3)
                vs3list.append(vis3)
        
            if row.find('td').get_text() == 'VIS 2':
                vis2 = row.find('td').next_element.next_element.get_text()
        
                if len(vis2) >= 10:
                    vis2l = vis2.split()
                    vis2 = float(vis2l[0]) * (10 ** float(vis2l[2][-1]))
        
                vis2 = float(vis2)
                vs2list.append(vis2)
        
            if row.find('td').get_text() == 'VIS 1':
                vis1 = row.find('td').next_element.next_element.get_text()
        
                if len(vis1) >= 10:
                    vis1l = vis1.split()
                    vis1 = float(vis1l[0]) * (10 ** float(vis1l[2][-1]))
        
                vis1 = float(vis1)
                vs1list.append(vis1)

            if row.find('td').get_text() == 'VIS ND1':
                nd1 = row.find('td').next_element.next_element.get_text()
        
                if len(nd1) >=10:
                    ndl = nd1.split()
                    nd1 = float(ndl[0]) * (10 ** float(ndl[2][-1]))
        
                nd1 = float(nd1)
                nd1list.append(nd1)
        
            if row.find('td').get_text() == 'VIS BK-7':
                bk7 = row.find('td').next_element.next_element.get_text()
        
                if len(bk7) >= 10:
                    bkl = bk7.split()
                    bk7 = float(bkl[0]) * (10 ** float(bkl[2][-1]))
        
                bk7 = float(bk7)
                bk7list.append(bk7)

    # To convert ra_deg and dec_deg to ra_hms and dec_dms.
    coord = SkyCoord(list(zip(ra_deg, dec_deg)), frame = 'icrs', unit = 'deg')
    RAhms_DECdms = coord.to_string('hmsdms', sep = ':')
    ra_hms, dec_dms = zip(*[hmdm.split(' ') for hmdm in RAhms_DECdms])

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

    all_together = np.array([ra_hms, dec_dms, cm, ct,
                            spty, vs3list, vs2list,
                            vs1list, nd1list, bk7list]).T

    # To select safe filters.
    filter_dict = {0: 'VIS3', 
                   1: 'VIS2', 
                   2: 'VIS1', 
                   3: 'ND1', 
                   4: 'BK7'}
    i = 0
    safe_filters = []
    for Filter in [vs3list, vs2list, vs1list, nd1list, bk7list]:
        if sum(Filter > 4800) == 0: 
            safe_filters.append(filter_dict[i])
        i = i + 1

    # To check if objects within 10 arcsec are present.
    sep = [ca.separation(cb).arcsecond for ca in coord for cb in coord 
           if ca.ra.value != cb.ra.value]

    proximity_check = np.array(sep) < proximity
    too_close = sum(proximity_check) / 2


    # Showing back the user inputs to user!
    print('\nPayload: {}, Coordinates: {}\n'.format(instrument, RADEC))

    #The usual mambo-jambo.
    print("\n\nTable of results")
    print("#########################\n")

    print("ra_hms\tdec_dms\tmag\tB-V\tSpecType\tVIS3\tVIS2\tVIS1\tND1\tBK7\n")
    for j in range(len(cm)):
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(
              ra_hms[j], dec_dms[j], cm[j], ct[j], spty[j], 
              vs3list[j], vs2list[j], vs1list[j], nd1list[j], 
              bk7list[j]))

    print("\n\nSafe filters: {}\n".format(safe_filters))

    if too_close > 0:
        print('\nWARNING! there exists {} pair of bright stars which are closer than\
              \n{} arcseconds!'.format(too_close, proximity))

    file_name = instrument + '_RA_' + RA_user + '_DEC_' + DEC_user + '.txt'
    with open(file_name,'w') as fr:
        fr.write("#ra_hms\tdec_dms\tmag\tB-V\tSpecType\tVIS3\tVIS2\tVIS1\tND1\tBK7\n")
        for j in range(len(cm)):
            fr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(
                     ra_hms[j], dec_dms[j], cm[j], ct[j], spty[j], 
                     vs3list[j], vs2list[j], vs1list[j], nd1list[j], 
                     bk7list[j]))  

        fr.write("\n\nSafe filters: {}\n".format(safe_filters))
        if too_close > 0:
            fr.write('\nWARNING! there exists {} pair of stars which are closer than\
                    \n{} arcseconds!'.format(too_close, proximity))

    print('\nDone\n')





