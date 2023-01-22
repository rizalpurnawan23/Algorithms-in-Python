# BISMILLAHIRRAHMANIRRAHIM-------------------------------------------
"""
Date        : 22/01/2023
Subject     : Email Extractor Algorithm
Authors     : Rizal Purnawan

This module contains an algorithm for extracting messages from
an email. The main purpose of the creation of this algorithm is to
ease document controling works in my office, so we can monitor
incoming and sent emails slightly authomatically. The algorithm uses
`imaplib` and `email` libraries in Python. This algorithm is still
in development and still needs a lot of improvements.

Rizal Purnawan
"""

def email_reader(
        the_email,
        email_server= "outlook.office365.com",
        port= 993,
        e_select= "inbox",
        e_search_par = "ALL",
        save= [None, False]):

    # Libraries:
    import imaplib
    import email
    import getpass
    import pandas as pd

    # Processing the email:
    mail = imaplib.IMAP4_SSL(email_server, port= port)
    # the_email = getpass.getpass("Email: ")
    the_email = the_email
    password = getpass.getpass("Password: ")
    mail.login(the_email, password)
    mail.select(e_select)
    _, message_num = mail.search(None, e_search_par)
    e_date = list()
    e_from = list()
    e_to = list()
    e_cc = list()
    e_bcc = list()
    e_subject = list()

    for mn in message_num[0].split():
        _, data = mail.fetch(mn, "(RFC822)")
        message = email.message_from_bytes(data[0][1])
        e_date.append(message.get('Date'))
        e_from.append(message.get('From'))
        e_to.append(message.get('To'))
        e_cc.append(message.get('CC'))
        e_bcc.append(message.get('BCC'))
        e_subject.append(message.get('Subject'))

    mail.close()

    # Creating a dataframe from the extracted emails:
    email_recap_dict = {
        "Date": e_date,
        "From": e_from,
        "To": e_to,
        "CC": e_cc,
        "BCC": e_bcc,
        "Subject": e_subject
    }
    email_recap_df = pd.DataFrame(email_recap_dict)
    if save[0] == None and save[1] == False:
        return email_recap_df
    elif type(save[0]) == str and save[1] == True:
        xl_name = save[0]
        if ".xlsx" in xl_name:
            pass
        else:
            xl_name = xl_name + ".xlsx"
        email_recap_df.to_excel(xl_name)
    else:
        print("ERROR: Improper 'save' parameters!")
        raise ValueError