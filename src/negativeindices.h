vector<OneBased> one_based_negative_indices(const SignConfiguration& v) {
    vector<OneBased> result;
    for (int i=1;i<=v.size();++i)
        if (v[i-1]<0) result.push_back(i);
    return result;
}
