function soln = interpolater(known1A, known2A, known1B, known2B, givenB)

    soln = known1A - ((known1A - known2A)*(known1B - givenB)/(known1B - known2B));

end