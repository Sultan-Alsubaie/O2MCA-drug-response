# data

| File | Content |
|---|---|
| `data.rds` | Reference results and annotation tables read by the scripts (32 files) |
| `inputs.rds.part00`–`part06` | One compressed bundle, split into parts below 25 MB, holding the four original CSV files and the two external fitting inputs |

Each script restores these files automatically the first time it runs. The restored files are byte-identical to the originals and are checked against the SHA-256 values in `config/`.
