//! The layout claim, asserted rather than assumed.
//!
//! `Alignment` is held O(n^2) live during an all-against-all pass, so a `Span`
//! that is wider than the two `u64` it replaces would be a real cost there.

use aln_coord::Span;
use std::mem::size_of;

#[test]
fn span_is_no_wider_than_the_two_integers_it_replaces() {
    assert_eq!(size_of::<Span>(), 2 * size_of::<u64>());
}

#[test]
fn one_optional_span_is_narrower_than_two_optional_integers() {
    // u64 has no niche, so a discriminant is unavoidable either way.  Span
    // pays for one instead of two.
    assert!(size_of::<Option<Span>>() < 2 * size_of::<Option<u64>>());
}
