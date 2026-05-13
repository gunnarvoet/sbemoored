# History

## Unreleased

### Breaking changes

- `sbemoored.sbe56.read_csv` now stores the human-readable short serial
  number in `attrs["SN"]` (e.g. `418`, `6413`) instead of the device-reported
  encoded value (e.g. `5_600_418`, `5_606_413`). The original device-form
  serial is preserved in the new `attrs["device_sn"]` for traceability.
  Code that previously joined `ds.SN` against a human-readable sensor
  database had to subtract `5_600_000` manually; that workaround is no
  longer needed. Existing NetCDFs written by older versions still carry
  the encoded value and either need to be regenerated or rewritten in
  place (`attrs["SN"] = old_SN - 5_600_000`, plus a new `attrs["device_sn"]
  = old_SN`).
