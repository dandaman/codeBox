## googleChromeCSV2pass.py

Little helper script to import Google Chrome Password.csv to Linux tool pass
Adapted from existing python tool to import Keepass2 csv

### Instructions
In your Google Chrome (or Chromium, whichever you use), type the following in your URL bar: chrome://flags/#password-import-export, and then enable the feature.
Restart your browser.
Go to your passwords chrome://settings/passwords (you may have to wait a little while for your passwords to sync), then scroll down to below your passwords and you'll see two new buttons, Import & Export.
Click Export, making sure you select the correct format (CSV).

`python3 ./googleChromeCSV2pass.py Google.csv [optional_folder_name|default=chromePW]`