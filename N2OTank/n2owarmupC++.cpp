// No Debugging has been done yet, but this is the C++ ver. 
                  }, 0.0, used);
            }
            t += used;
            T = next;
            if (result.full) {
                result.points.emplace_back(t, T);
                if (std::abs(t - 900.0) < 1e-7) result.T15 = T;
                return result;
            }
            dt *= factor;
        }
        result.points.emplace_back(t, T);
        if (sample == 90) result.T15 = T;
    }
    return result;
}

std::string fmt(double value, int precision = 1) {
    if (!std::isfinite(value)) return "n/a";
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}
struct Curve { std::string name; Result result; };
void save_outputs(const std::vector<Curve>& curves, double T_full) {
    std::ofstream csv("n2o_tank_warmup.csv"), svg("n2o_tank_warmup.svg");
    if (!csv || !svg) throw std::runtime_error("Cannot open output files.");
    csv << "case,time_min,temperature_K\n" << std::setprecision(12);
    double low = T0, high = std::max(T_full, T_goal);
    for (const auto& c : curves) for (const auto& p : c.result.points) {
        csv << '"' << c.name << "\"," << p.first / 60.0 << ',' << p.second << '\n';
        low = std::min(low, p.second); high = std::max(high, p.second);
    }
    low = std::floor(low - 1.0); high = std::ceil(high + 2.0);
    const auto x = [](double minutes) { return 75.0 + 620.0 * minutes / 30.0; };
    const auto y = [&](double T) { return 350.0 - 290.0 * (T - low) / (high - low); };
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"760\" height=\"470\" viewBox=\"0 0 760 470\">"
        << "<rect width=\"760\" height=\"470\" fill=\"white\"/>"
        << "<g font-family=\"sans-serif\" font-size=\"12\">"
        << "<text x=\"380\" y=\"25\" text-anchor=\"middle\">Exposed 6 in x 0.25 Al tank, IREC sun, 5.2 m/s wind</text>";
    for (int m = 0; m <= 30; m += 5)
        svg << "<path d=\"M" << x(m) << " 60V350\" stroke=\"#ddd\"/>"
            << "<text x=\"" << x(m) << "\" y=\"370\" text-anchor=\"middle\">" << m << "</text>";
    for (double T = low; T <= high; T += 2.0)
        svg << "<path d=\"M75 " << y(T) << "H695\" stroke=\"#ddd\"/>"
            << "<text x=\"65\" y=\"" << y(T) + 4 << "\" text-anchor=\"end\">" << T << "</text>";
    svg << "<path d=\"M75 60V350H695\" fill=\"none\" stroke=\"black\"/>"
        << "<path d=\"M75 " << y(T_goal) << "H695\" stroke=\"black\"/>"
        << "<path d=\"M" << x(15) << " 60V350\" stroke=\"gray\" stroke-dasharray=\"3 3\"/>"
        << "<path d=\"M75 " << y(T_full) << "H695\" stroke=\"red\" stroke-dasharray=\"3 3\"/>"
        << "<text x=\"80\" y=\"" << y(T_full) - 5 << "\" fill=\"red\">liquid-full</text>"
        << "<text x=\"385\" y=\"395\" text-anchor=\"middle\">time after fill [min]</text>"
        << "<text transform=\"translate(20 205) rotate(-90)\" text-anchor=\"middle\">N2O temperature [K]</text>";
    const std::vector<std::string> colors = {"#1675b6", "#e87916"};
    for (std::size_t i = 0; i < curves.size(); ++i) {
        const auto& c = curves[i];
        svg << "<polyline fill=\"none\" stroke=\"" << colors[i % colors.size()]
            << "\" stroke-width=\"2\" points=\"";
        for (const auto& p : c.result.points) svg << x(p.first / 60.0) << ',' << y(p.second) << ' ';
        svg << "\"/><text x=\"90\" y=\"" << 420 + 20 * i << "\" fill=\""
            << colors[i % colors.size()] << "\">" << c.name << "</text>";
    }
    svg << "</g></svg>\n";
    csv.flush(); svg.flush();
    if (!csv || !svg) throw std::runtime_error("Failed to write output files.");
}

int main() {
    try {
        const Model model;
        const std::vector<std::pair<std::string, double>> T_AIR = {
            {"307.15 K air", 307.15}, {"293.65 K air", 293.65}};
        struct Finish { std::string name; double alpha, eps; };
        const std::vector<Finish> FINISH = {{"bare Al", 0.35, 0.09}};
        std::cout << "per meter: m_Al = " << fmt(m_Al, 3)
            << " kg/m, m_N2O = " << fmt(model.m_N2O, 3)
            << " kg/m (rho_bar " << fmt(model.rho_bar) << " kg/m3), c_eff = "
            << fmt((model.C_closed_form - m_Al * c_Al) / model.m_N2O, 0)
            << " J/kg-K, C' = " << fmt(model.C_closed_form, 0)
            << " J/m-K, liquid-full at " << fmt(model.T_full) << " K\n\n";
        const std::vector<std::string> cols = {"Tank, Air Temp", "h_c", "h_r", "T_sa",
            "UA'", "tau", "t288 est", "t288 ODE", "T15", "T30", "T0 for 288@15"};
        std::vector<int> widths;
        for (const auto& c : cols) widths.push_back(std::max(7, static_cast<int>(c.size())));
        widths[0] = 22; widths[9] = 18;
        const auto row = [&](const std::vector<std::string>& values) {
            for (std::size_t i = 0; i < values.size(); ++i)
                std::cout << (i == 0 ? std::left : std::right) << std::setw(widths[i])
                          << values[i] << (i + 1 == values.size() ? "\n" : "  ");
        };
        row(cols);
        std::vector<Curve> curves;
        for (const auto& air : T_AIR) for (const auto& fin : FINISH) {
            const auto a = air_side(0.5 * (T0 + T_goal), air.second, fin.alpha, fin.eps);
            const double tau = model.C_closed_form / a.UA;
            double t_est = std::numeric_limits<double>::quiet_NaN();
            if (a.T_sa > T_goal)
                t_est = tau * std::log((a.T_sa - T0) / (a.T_sa - T_goal)) / 60.0;
            const double T0_req = a.T_sa - (a.T_sa - T_goal) * std::exp(900.0 / tau);
            const Result result = integrate(model, air.second, fin.alpha, fin.eps);
            const double stop = result.points.back().first;
            const std::string name = fin.name + ", " + air.first;
            row({name, fmt(a.h_c), fmt(a.h_r), fmt(a.T_sa), fmt(a.UA, 2),
                 fmt(tau / 60.0), fmt(t_est), fmt(result.hit / 60.0), fmt(result.T15),
                 result.full ? "full @ " + fmt(stop / 60.0, 0) + " min" :
                               fmt(result.points.back().second), fmt(T0_req)});
            curves.push_back({name, result});
        }
        std::cout << "\nunits: h W/m2-K | T K | UA' W/m-K | tau, t288 min\n"
                  << "T0 for 288@15 uses the constant-capacity estimate, as in Python.\n"
                  << "n/a: target not reached, or requested time after liquid-full stop.\n";
        save_outputs(curves, model.T_full);
        std::cout << "saved: n2o_tank_warmup.csv, n2o_tank_warmup.svg\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}
