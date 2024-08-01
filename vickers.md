vickers - a variantkey implementation in Rust
==

A Rust implementation of [variantkey](https://github.com/Genomicsplc/variantkey).

This is a rather `panic!`y implementation and not necessarily fit for production purposes. It's mostly for my benefit
to try out Rust.

Once you have done `cargo install --path ./` you can run the CLI like so:

```shell
vickers --help
```

Convert a variant to a key:

```shell
vickers vk 10-123456-T-G
# 5764872642925428736
```
output in hexadecimal:

```shell
vickers vk -x 10-123456-T-G
# 5000f12008f00000
```

And a key back to a variant:

```shell
vickers kv 5764872642925428736
# 10-123456-T-G
```
or from hexadecimal:

```shell
vickers kv -x 5000f12008f00000
# 10-123456-T-G
```