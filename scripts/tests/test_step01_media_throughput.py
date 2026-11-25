#!/usr/bin/env python3
"""Unit tests for Step 01 media throughput: wave rescue, host interleaving, circuit breaker, QC cache, and GBIF dedup."""

from __future__ import annotations

import asyncio
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
from PIL import Image

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.core.harvester as _harv
from scripts.core.harvester import VoucherHarvester, download_all_voucher_images, download_single_image
from scripts.core.step01_media_qc_cache import (
    audit_voucher_images_parallel,
    encode_and_write_image_bytes,
    interleave_records_by_host,
)


def _make_jpeg_bytes(width: int = 3200, height: int = 3000) -> bytes:
    buf = io.BytesIO()
    Image.new("RGB", (width, height), color=(115, 130, 108)).save(buf, format="JPEG", quality=92)
    return buf.getvalue()


class TestStep01MediaThroughputOptimizations(unittest.TestCase):
    """Verify Step 01 wave-based concurrent rescue, host scheduling, dead-host fast-fail, QC cache, and GBIF dedup."""

    def test_interleave_records_by_host_eliminates_host_clustering(self) -> None:
        """Round-robin host interleaving spreads consecutive tasks across distinct servers."""
        clustered = [
            ("https://host-a.org/1.jpg", Path("a1.webp")),
            ("https://host-a.org/2.jpg", Path("a2.webp")),
            ("https://host-a.org/3.jpg", Path("a3.webp")),
            ("https://host-b.org/1.jpg", Path("b1.webp")),
            ("https://host-b.org/2.jpg", Path("b2.webp")),
            ("https://host-c.org/1.jpg", Path("c1.webp")),
        ]
        interleaved = interleave_records_by_host(clustered)
        hosts = [u.split("/")[2] for u, _ in interleaved]
        self.assertEqual(hosts[:3], ["host-a.org", "host-b.org", "host-c.org"])
        self.assertEqual(hosts[3:5], ["host-a.org", "host-b.org"])
        self.assertEqual(hosts[5], "host-a.org")

    def test_encode_and_write_image_bytes_webp_and_passthrough(self) -> None:
        """JPEG payloads convert cleanly to WebP (method=2), while native WebP bytes pass through untouched."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            jpg_bytes = _make_jpeg_bytes(400, 300)
            dest_webp = tmp_p / "NCU_TEST_01.webp"
            self.assertTrue(encode_and_write_image_bytes(jpg_bytes, dest_webp, webp_quality=90, image_mod=Image))
            self.assertTrue(dest_webp.is_file())
            written = dest_webp.read_bytes()
            self.assertEqual(written[:4], b"RIFF")
            self.assertEqual(written[8:12], b"WEBP")

            dest_pass = tmp_p / "NCU_TEST_02.webp"
            self.assertTrue(encode_and_write_image_bytes(written, dest_pass, webp_quality=90, image_mod=Image))
            self.assertEqual(dest_pass.read_bytes(), written)

    def test_download_fast_fail_permanent_errors_and_dead_host_circuit_breaker(self) -> None:
        """HTTP 403/404 and HTML 200 fail on attempt 1; dead host trips circuit breaker after 6 consecutive network errors."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            calls_by_url: dict[str, int] = {}

            class _Resp:
                def __init__(self, status: int, content_type: str, body: bytes) -> None:
                    self.status = status
                    self.headers = {"Content-Type": content_type}
                    self._body = body

                async def read(self) -> bytes:
                    return self._body

                async def __aenter__(self) -> "_Resp":
                    return self

                async def __aexit__(self, *_a) -> None:
                    return None

            class _TimeoutCtx:
                async def __aenter__(self):
                    raise asyncio.TimeoutError("connect timeout")

                async def __aexit__(self, *_a) -> None:
                    return None

            valid_jpg = b"\xff\xd8\xff" + b"X" * 2048

            class _FakeSession:
                def get(self, url: str, timeout=None):
                    calls_by_url[url] = calls_by_url.get(url, 0) + 1
                    if "forbidden" in url:
                        return _Resp(403, "text/html", b"forbidden")
                    if "html200" in url:
                        return _Resp(200, "text/html; charset=utf-8", b"<html>not an image</html>" * 100)
                    if "dead.herbarium.test" in url:
                        return _TimeoutCtx()
                    return _Resp(200, "image/jpeg", valid_jpg)

            sem = asyncio.Semaphore(4)
            session = _FakeSession()

            # 1. Permanent 403 and HTML 200 fail immediately on attempt 1 (zero retries)
            ok_403, _ = asyncio.run(download_single_image(session, "https://ok.test/forbidden.jpg", tmp_p / "a.jpg", sem))
            ok_html, _ = asyncio.run(download_single_image(session, "https://ok.test/html200.jpg", tmp_p / "b.jpg", sem))
            self.assertFalse(ok_403)
            self.assertFalse(ok_html)
            self.assertEqual(calls_by_url["https://ok.test/forbidden.jpg"], 1)
            self.assertEqual(calls_by_url["https://ok.test/html200.jpg"], 1)

            # 2. Dead host circuit breaker trips after 6 failed URLs and fast-fails URLs #7..#12 with 0 network calls
            host_health: dict[str, dict[str, int]] = {}
            with patch("scripts.core.step01_harvester_media.asyncio.sleep", return_value=None):
                for idx in range(12):
                    u_dead = f"https://dead.herbarium.test/img_{idx}.jpg"
                    ok_d, _ = asyncio.run(
                        download_single_image(
                            session, u_dead, tmp_p / f"d_{idx}.jpg", sem,
                            max_retries=3, host_health=host_health,
                        )
                    )
                    self.assertFalse(ok_d)
                # Healthy host still succeeds cleanly in the same health table
                ok_good, _ = asyncio.run(
                    download_single_image(
                        session, "https://alive.herbarium.test/good.jpg", tmp_p / "good.jpg", sem,
                        host_health=host_health,
                    )
                )
            self.assertTrue(ok_good)
            self.assertEqual(host_health["dead.herbarium.test"]["tripped"], 1)
            for idx in range(6, 12):
                self.assertNotIn(f"https://dead.herbarium.test/img_{idx}.jpg", calls_by_url)

    def test_wave_based_concurrent_secondary_duplicate_rescue(self) -> None:
        """Secondary duplicate rescue batches unresolved exsiccatae into concurrent waves (never concurrency=1)."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "data" / "raw_vouchers"
            raw_dir.mkdir(parents=True, exist_ok=True)
            good_bytes = _make_jpeg_bytes(3200, 3000)

            harv = VoucherHarvester(
                taxa=["Packera anonyma"],
                workspace_dir=tmp_p,
                raw_dir=raw_dir,
                min_megapixels=4.0,
                min_file_size_kb=5.0,
                check_sharpness=False,
                concurrency=16,
                expand_synonyms=False,
            )
            df_prim = pd.DataFrame([
                {
                    "catalogNumber": "PRIM_1", "institutionCode": "NCU", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "101", "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold", "exsiccatae_key": "barkley_101_1980-05-10",
                    "image_path": "data/raw_vouchers/PRIM_1.webp", "_image_url": "https://ex.test/prim_1_bad.jpg",
                },
                {
                    "catalogNumber": "PRIM_2", "institutionCode": "NCU", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "102", "eventDate": "1980-05-11",
                    "determiner_tier": "Tier_1_Gold", "exsiccatae_key": "barkley_102_1980-05-11",
                    "image_path": "data/raw_vouchers/PRIM_2.webp", "_image_url": "https://ex.test/prim_2_bad.jpg",
                },
            ])
            harv._secondary_duplicate_pool = pd.DataFrame([
                # ek1 rank-1 Gold succeeds in Wave 1; ek1 rank-2 Bronze must NEVER be downloaded
                {
                    "catalogNumber": "SEC_1_GOLD", "institutionCode": "MO", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "101", "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold", "exsiccatae_key": "barkley_101_1980-05-10",
                    "image_path": "data/raw_vouchers/SEC_1_GOLD.webp", "_image_url": "https://ex.test/sec_1_gold_ok.jpg",
                },
                {
                    "catalogNumber": "SEC_1_BRZ", "institutionCode": "GA", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "101", "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_3_Bronze", "exsiccatae_key": "barkley_101_1980-05-10",
                    "image_path": "data/raw_vouchers/SEC_1_BRZ.webp", "_image_url": "https://ex.test/sec_1_brz_unused.jpg",
                },
                # ek2 rank-1 Silver fails in Wave 1; ek2 rank-2 Bronze succeeds in Wave 2
                {
                    "catalogNumber": "SEC_2_SLV", "institutionCode": "US", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "102", "eventDate": "1980-05-11",
                    "determiner_tier": "Tier_2_Silver", "exsiccatae_key": "barkley_102_1980-05-11",
                    "image_path": "data/raw_vouchers/SEC_2_SLV.webp", "_image_url": "https://ex.test/sec_2_slv_bad.jpg",
                },
                {
                    "catalogNumber": "SEC_2_BRZ", "institutionCode": "NY", "scientificName": "Packera anonyma",
                    "recordedBy": "Barkley", "recordNumber": "102", "eventDate": "1980-05-11",
                    "determiner_tier": "Tier_3_Bronze", "exsiccatae_key": "barkley_102_1980-05-11",
                    "image_path": "data/raw_vouchers/SEC_2_BRZ.webp", "_image_url": "https://ex.test/sec_2_brz_ok.jpg",
                },
            ])

            recorded_calls: list[tuple[int, list[str]]] = []

            async def _fake_batch_dl(records_to_download, concurrency_limit=15, **_kw):
                urls = [str(u) for u, _ in records_to_download]
                recorded_calls.append((int(concurrency_limit), urls))
                succ = 0
                for u, d in records_to_download:
                    if "_ok" in str(u):
                        encode_and_write_image_bytes(good_bytes, Path(d), webp_quality=92, image_mod=Image)
                        succ += 1
                return {"success": succ, "skipped": 0, "failed": len(records_to_download) - succ}

            with patch("scripts.core.harvester.download_all_voucher_images", side_effect=_fake_batch_dl):
                df_out, stats_out = harv.download_and_validate_media(df_prim)

            # Exactly 3 batch calls: Primary + Wave 1 (both rank-1 candidates together) + Wave 2 (only ek2 rank-2)
            self.assertEqual(len(recorded_calls), 3)
            self.assertEqual([c for c, _ in recorded_calls], [16, 16, 16])
            self.assertEqual(set(recorded_calls[1][1]), {"https://ex.test/sec_1_gold_ok.jpg", "https://ex.test/sec_2_slv_bad.jpg"})
            self.assertEqual(recorded_calls[2][1], ["https://ex.test/sec_2_brz_ok.jpg"])
            self.assertEqual(set(df_out["catalogNumber"]), {"SEC_1_GOLD", "SEC_2_BRZ"})
            self.assertEqual(stats_out["quality_rejected"], 2)

    def test_parallel_qc_audit_cache_hit_threshold_reeval_and_mtime_invalidation(self) -> None:
        """Cached QC raw metrics avoid re-decoding unchanged files, re-evaluate new thresholds, and invalidate on change."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            img_dir = tmp_p / "raw_vouchers"
            img_dir.mkdir(parents=True, exist_ok=True)
            cache_p = tmp_p / ".voucher_image_qc_cache.json"

            paths = []
            for idx in range(6):
                p = img_dir / f"V_{idx}.webp"
                encode_and_write_image_bytes(_make_jpeg_bytes(3200, 3000), p, webp_quality=90, image_mod=Image)
                paths.append(p)

            res_cold = audit_voucher_images_parallel(
                paths, min_megapixels=8.0, min_file_size_kb=5.0, check_sharpness=False,
                max_workers=4, cache_path=cache_p,
            )
            self.assertTrue(all(ok for ok, _ in res_cold))
            self.assertTrue(cache_p.is_file())

            # Warm run with stricter MP threshold (15.0 > 9.6 MP) uses cache and rejects cleanly without PIL open
            res_warm_strict = audit_voucher_images_parallel(
                paths, min_megapixels=15.0, min_file_size_kb=5.0, check_sharpness=False,
                max_workers=4, cache_path=cache_p,
            )
            self.assertTrue(all((not ok and m["reason"] == "low_resolution") for ok, m in res_warm_strict))

            # Modify file #0 to low-res (change size/mtime) -> cache invalidates and detects new size
            encode_and_write_image_bytes(_make_jpeg_bytes(400, 400), paths[0], webp_quality=90, image_mod=Image)
            st0 = paths[0].stat()
            os.utime(paths[0], ns=(st0.st_atime_ns, st0.st_mtime_ns + 5_000_000))
            res_after_mod = audit_voucher_images_parallel(
                paths, min_megapixels=8.0, min_file_size_kb=1.0, check_sharpness=False,
                max_workers=4, cache_path=cache_p,
            )
            self.assertFalse(res_after_mod[0][0])
            self.assertTrue(all(ok for ok, _ in res_after_mod[1:]))

    def test_gbif_duplicate_species_key_synonym_short_circuits_pagination(self) -> None:
        """Live pagination stops at page 1 and skips hasCoordinate=False when a synonym maps to an exhausted speciesKey."""
        page0 = [
            {
                "key": 80000 + i, "speciesKey": 3109240, "scientificName": "Packera anonyma",
                "stateProvince": "North Carolina", "decimalLatitude": 35.8 + i * 0.001, "decimalLongitude": -79.0,
                "coordinateUncertaintyInMeters": 50.0, "year": 2020, "month": 5, "day": 10,
                "recordedBy": "Collector", "recordNumber": str(i), "identifiedBy": "T.M. Barkley",
                "media": [{"identifier": f"https://ex.test/a_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": f"NCU_A_{i}", "institutionCode": "NCU",
            }
            for i in range(4)
        ]
        # Add one rejected record (missing year) so visited_country_gbif_keys (not just retained keys) is tested
        page0.append({
            "key": 89999, "speciesKey": 3109240, "scientificName": "Packera anonyma",
            "stateProvince": "North Carolina", "decimalLatitude": 35.85, "decimalLongitude": -79.05,
            "media": [{"identifier": "https://ex.test/no_year.jpg", "type": "StillImage", "format": "image/jpeg"}],
            "catalogNumber": "NCU_NO_YEAR", "institutionCode": "NCU",
        })
        page1 = [
            {
                "key": 81000 + i, "speciesKey": 3109240, "scientificName": "Packera anonyma",
                "stateProvince": "North Carolina", "decimalLatitude": 35.9 + i * 0.001, "decimalLongitude": -79.1,
                "coordinateUncertaintyInMeters": 50.0, "year": 2021, "month": 5, "day": 12,
                "recordedBy": "Collector", "recordNumber": str(100 + i), "identifiedBy": "T.M. Barkley",
                "media": [{"identifier": f"https://ex.test/b_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": f"NCU_B_{i}", "institutionCode": "NCU",
            }
            for i in range(2)
        ]

        seen_calls: list[dict] = []

        class _LiveLikeOcc:
            @staticmethod
            def search(**kw):
                seen_calls.append(dict(kw))
                if kw.get("typeStatus") == "*":
                    return {"results": [], "count": 0}
                if kw.get("hasCoordinate") is False:
                    return {"results": [], "count": 0}
                if kw.get("country") == "CA":
                    return {"results": [], "count": 0}
                off = int(kw.get("offset", 0))
                if off == 0:
                    return {"results": page0, "count": len(page0) + len(page1)}
                return {"results": page1, "count": len(page0) + len(page1)}

        orig_occ = _harv.occ
        _harv.occ = _LiveLikeOcc()
        try:
            with patch("scripts.core.harvester.time.sleep", return_value=None), \
                 patch("scripts.core.harvester.filter_vouchers_by_usda_county_range", side_effect=lambda df, **_: (df, {})):
                harv = VoucherHarvester(
                    taxa=["Packera anonyma", "Senecio smallii"],
                    countries=["US", "CA"],
                    max_records_per_taxon=50,
                    expand_synonyms=False,
                    include_ungeoreferenced_state_vouchers=True,
                )
                df_h = harv.harvest()
        finally:
            _harv.occ = orig_occ

        self.assertEqual(len(df_h), 6)
        syn_calls = [c for c in seen_calls if c.get("scientificName") == "Senecio smallii"]
        syn_us_offsets = [c.get("offset") for c in syn_calls if c.get("country") == "US" and "typeStatus" not in c and c.get("hasCoordinate") is not False]
        syn_ungeo_calls = [c for c in syn_calls if c.get("hasCoordinate") is False]
        self.assertEqual(syn_us_offsets, [0], "Duplicate speciesKey synonym must stop at offset=0 without re-paginating")
        self.assertEqual(syn_ungeo_calls, [], "Duplicate speciesKey synonym must skip redundant hasCoordinate=False pass")


if __name__ == "__main__":
    unittest.main()

