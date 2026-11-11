"""What a store holds beside a product's bytes: the names the orchestration's reconciler reads.

The toolkit imports nothing from ``auto_research`` (CLAUDE.md §1), so the store's product contract is
written here once, and ``auto_research``'s test suite checks that these values equal the reconciler's
(``auto_research.store.reconcile``: ``META_*``, ``RECIPE_FILE``, ``RecipeFile``/``RecipeObject``).

* Every object a product stores carries five immutable metadata entries (``x-amz-meta-<name>``):
  :data:`META_NAMES`.
* A product's folder ends with :data:`RECIPE_FILE`, written last: the recipe as claimed, its hash, the
  claim generation, and every other object of the product (:func:`recipe_file`).
* A recipe holds exactly :data:`RECIPE_KEYS`; each input is a record with :data:`RECORD_KEYS`.
"""

from __future__ import annotations

import re

META_RECIPE_HASH = "recipe-hash"
META_SHA256 = "sha256"
META_PRODUCER = "producer"
META_TOOLKIT_COMMIT = "toolkit-commit"
META_ENV_FINGERPRINT = "env-fingerprint"
#: The metadata every stored object carries, as boto3 names them (``Metadata`` keys).
META_NAMES = (META_RECIPE_HASH, META_SHA256, META_PRODUCER, META_TOOLKIT_COMMIT, META_ENV_FINGERPRINT)

#: The completeness marker of a product's folder.
RECIPE_FILE = "recipe.json"
#: The keys of a recipe (catalog.claim refuses any other set).
RECIPE_KEYS = ("kind", "subject", "tool", "algorithm_version", "settings", "inputs")
#: The keys of a record: an object version in a named store.
RECORD_KEYS = ("store", "bucket", "key", "version_id")
#: The keys of :data:`RECIPE_FILE`, and of each object it lists.
RECIPE_FILE_KEYS = ("recipe_hash", "generation", "recipe", "objects")
RECIPE_OBJECT_KEYS = ("key", "version_id", "sha256", "size")

_HEX64 = re.compile(r"^[0-9a-f]{64}$")


class ContractError(ValueError):
    """A recipe file or record that does not have the contract's shape."""


def recipe_file(recipe_hash: str, generation: int, recipe: dict, objects: list[dict]) -> dict:
    """The body of :data:`RECIPE_FILE`, checked: ``objects`` lists every other object of the product
    as ``{key, version_id, sha256, size}``, at least one."""
    body = {"recipe_hash": recipe_hash, "generation": generation, "recipe": recipe, "objects": objects}
    check_recipe_file(body)
    return body


def check_recipe_file(body) -> dict:
    """``body`` when it has the shape of :data:`RECIPE_FILE`; raises :class:`ContractError` naming what
    differs."""
    if not isinstance(body, dict) or sorted(body) != sorted(RECIPE_FILE_KEYS):
        raise ContractError(f"a {RECIPE_FILE} holds exactly {', '.join(RECIPE_FILE_KEYS)}")
    if not isinstance(body["recipe_hash"], str) or not _HEX64.match(body["recipe_hash"]):
        raise ContractError(f"{RECIPE_FILE}: recipe_hash is 64 lower-case hex")
    if type(body["generation"]) is not int or body["generation"] < 1:
        raise ContractError(f"{RECIPE_FILE}: generation is a whole number from 1")
    if not isinstance(body["recipe"], dict) or sorted(body["recipe"]) != sorted(RECIPE_KEYS):
        raise ContractError(f"{RECIPE_FILE}: recipe holds exactly {', '.join(RECIPE_KEYS)}")
    objects = body["objects"]
    if not isinstance(objects, list) or not objects:
        raise ContractError(f"{RECIPE_FILE}: objects lists at least one object")
    for item in objects:
        if not isinstance(item, dict) or sorted(item) != sorted(RECIPE_OBJECT_KEYS):
            raise ContractError(f"{RECIPE_FILE}: each object is {{{', '.join(RECIPE_OBJECT_KEYS)}}}")
        if not (isinstance(item["key"], str) and item["key"] and isinstance(item["version_id"], str)
                and item["version_id"] and isinstance(item["sha256"], str) and _HEX64.match(item["sha256"])
                and type(item["size"]) is int and item["size"] >= 0):
            raise ContractError(f"{RECIPE_FILE}: object {item!r} has a key, a version_id, a 64-hex "
                                "sha256 and a size in bytes")
    return body
